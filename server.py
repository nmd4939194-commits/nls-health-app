import os
import re
import json
import struct
import datetime
import requests
from typing import Optional, Dict, Any
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="NLS Health Sync Engine", version="4.5.0")

# 允許跨域請求，讓 Netlify 的前端網頁可以呼叫 API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DELPHI_BASE_DATE = datetime.datetime(1899, 12, 30)
RECORD_MAGIC_PREFIX = b'\xc8\xd1\xd1\xcb\xc5\xc4\xce'  # ИССЛЕДО (器官區塊開頭)

ORGAN_MAP_ZH = {
    "human19.bmp": "腹膜後腔器官",
    "5125_l.bmp": "胃前壁組織微觀結構",
    "version5.bmp": "小腸黏膜與胃壁微觀結構",
    "r30.bmp": "小腸壁組織微觀結構",
    "5124_l.bmp": "食道、胃、十二指腸綜合結構",
    "the liver, viewed from above.bmp": "肝臟解剖(上方)",
    "right lung (medial view).bmp": "右肺與內側氣道",
    "left lung (medial view).bmp": "左肺與內側氣道",
    "the kidney (left).bmp": "左腎縱切面",
    "the kidney (right).bmp": "右腎縱切面",
    "r47.bmp": "攝護腺組織",
    "cor1.bmp": "心臟前壁血管與冠狀動脈",
    "r16a.bmp": "動脈血管壁微循環",
    "5522.bmp": "腦底動脈環(威利氏環)",
    "66_65_r.bmp": "右側眼球結構",
    "66_65_l.bmp": "左側眼球結構",
    "the fundus of the right eye.bmp": "右眼底微循環",
    "the fundus of the right eye2.bmp": "右眼底黃斑部",
    "r1901.bmp": "脊椎椎體正中矢狀切面",
    "r03.bmp": "紅骨髓造血機能",
    "xromc.bmp": "C組染色體",
    "xrome.bmp": "E組染色體",
    "r3371.bmp": "子宮縱切面",
    "r22.bmp": "乳腺組織"
}

def extract_patient_info(raw_bytes: bytes) -> Dict[str, str]:
    head_slice = raw_bytes[:2000]
    head_latin = head_slice.decode('latin1', errors='ignore')
    
    patient_name = "受檢者"
    patient_dob = "未知"

    name_match = re.search(r'(?:tw|TW)([\xa1-\xfe][\x40-\x7e\xa1-\xfe]+)', head_latin)
    if name_match:
        try:
            name_bytes = name_match.group(1).encode('latin1')
            patient_name = name_bytes.decode('big5', errors='ignore').strip()
        except: pass
    else:
        for k in range(0, min(500, len(head_slice) - 6)):
            try:
                candidate = head_slice[k:k+8].decode('big5')
                if len(candidate) >= 2 and all('\u4e00' <= c <= '\u9fff' for c in candidate):
                    patient_name = candidate.strip()
                    break
            except: continue

    dob_match = re.search(r'(\d{4})(\d{2})(\d{2})(?:abc)?', head_latin)
    if dob_match:
        patient_dob = f"{dob_match.group(1)}-{dob_match.group(2)}-{dob_match.group(3)}"

    return { "name": patient_name, "dob": patient_dob }

def extract_organ_name(header_area: bytes, index: int) -> str:
    zh_name = ""
    # 1. 嘗試從雙破折號後解析 Big5 中文
    dash_pos = header_area.find(b'--')
    if dash_pos != -1:
        tail = header_area[dash_pos+2:dash_pos+120]
        match = re.search(rb'[A-Z]', tail)
        if match:
            zh_bytes = tail[:match.start()]
        else:
            zh_bytes = tail
            
        while zh_bytes and zh_bytes[-1] < 0x80 and chr(zh_bytes[-1]) not in '()（）-，、':
            zh_bytes = zh_bytes[:-1]
            
        try:
            parsed_zh = zh_bytes.decode('big5', errors='ignore').strip()
            if len(parsed_zh) >= 2: zh_name = parsed_zh
        except: pass

    # 2. 尋找 bmp 檔名並對應字典
    bmp_match = re.search(rb'([a-zA-Z0-9_\-\.\s\(\)]+\.[bB][mM][pP])', header_area)
    if bmp_match:
        bmp_file = bmp_match.group(1).decode('latin1', errors='ignore').strip().lower()
        bmp_file = re.sub(r'^_+', '', bmp_file)
        for key, val in ORGAN_MAP_ZH.items():
            if key in bmp_file or bmp_file in key:
                if not zh_name: zh_name = val
                break
                
    # 3. 嘗試提取英文名稱作為備案
    if not zh_name:
        en_candidates = re.findall(rb'[A-Z][A-Za-z0-9\s,\-\(\)\/\.]{4,}', header_area)
        for cand in en_candidates:
            cand_str = cand.decode('latin1', errors='ignore').strip()
            if not cand_str.lower().endswith('.bmp'):
                zh_name = cand_str
                break

    return zh_name if zh_name else f"檢測項目 #{index}"

def extract_entropy_counts(block_bytes: bytes):
    counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    total = 0
    
    # 策略一：尋找原機統計數據區塊 'ДАННЫЕ:'
    magic_data = b'\xc4\xc0\xcd\xcd\xdb\xc5\x3a'
    idx = block_bytes.find(magic_data)
    if idx != -1:
        search_area = block_bytes[idx + len(magic_data) : idx + len(magic_data) + 40]
        best_total = 0
        best_counts = None
        for offset in range(12): # 容錯偏移量
            try:
                ints = struct.unpack('<6i', search_area[offset:offset+24])
                # 驗證數值合理性
                if all(0 <= v < 50000 for v in ints):
                    cur_total = sum(ints)
                    if 5 <= cur_total < 100000 and cur_total > best_total:
                        best_total = cur_total
                        best_counts = {i+1: ints[i] for i in range(6)}
            except: pass
            
        if best_counts:
            return best_counts, best_total

    # 策略二：4相位滑動視窗容錯掃描 (應對無標頭器官)
    best_fallback_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    best_fallback_total = 0
    
    # 跳過標頭，由後段 3/4 開始精準捕捉 [值, 00, 00, 00]
    start_offset = len(block_bytes) // 4
    
    for align in range(4):
        cur_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        for off in range(start_offset + align, len(block_bytes) - 3, 4):
            val = block_bytes[off]
            if 1 <= val <= 6 and block_bytes[off+1] == 0 and block_bytes[off+2] == 0 and block_bytes[off+3] == 0:
                cur_counts[val] += 1
                
        cur_total = sum(cur_counts.values())
        if cur_total > best_fallback_total:
            best_fallback_total = cur_total
            best_fallback_counts = cur_counts
            
    return best_fallback_counts, best_fallback_total

def parse_single_record(block_bytes: bytes, index: int) -> Optional[Dict[str, Any]]:
    if len(block_bytes) < 64: return None

    # 解析時間
    record_time = "未知"
    for offset in range(8, min(120, len(block_bytes) - 8)):
        try:
            val = struct.unpack('<d', block_bytes[offset:offset+8])[0]
            if 40000.0 <= val <= 55000.0: 
                dt = DELPHI_BASE_DATE + datetime.timedelta(days=val)
                record_time = dt.strftime("%Y-%m-%d %H:%M")
                break
        except: continue

    zh_name = extract_organ_name(block_bytes[:800], index)
    entropy_counts, total_measured_points = extract_entropy_counts(block_bytes)
            
    weighted_sum = sum(k * v for k, v in entropy_counts.items())
    avg_entropy = round(weighted_sum / total_measured_points, 2) if total_measured_points > 0 else 2.5

    return {
        "index": index,
        "date_time": record_time,
        "zh_title": zh_name,
        "total_points": total_measured_points,
        "avg_entropy": avg_entropy,
        "entropy": {
            "level_1": entropy_counts[1], "level_2": entropy_counts[2], "level_3": entropy_counts[3],
            "level_4": entropy_counts[4], "level_5": entropy_counts[5], "level_6": entropy_counts[6]
        }
    }

@app.post("/api/upload-and-split") # 相容舊版名稱，保留原接口
@app.post("/api/upload-and-sync")
async def upload_and_sync(
    file: UploadFile = File(...),
    gas_url: str = Form("") 
):
    try:
        content = await file.read()
        patient = extract_patient_info(content)
        
        # 安全切割檔案，不漏接任何一個器官
        offsets = [m.start() for m in re.finditer(re.escape(RECORD_MAGIC_PREFIX), content)]
        offsets.append(len(content))
        
        records = []
        for i in range(len(offsets) - 1):
            chunk = content[offsets[i]:offsets[i+1]]
            item = parse_single_record(chunk, i + 1)
            # 過濾掉雜訊 (小於5點通常是圖檔雜訊)
            if item and item["total_points"] >= 5: 
                records.append(item)

        if not records:
            raise ValueError("無法解析出有效的器官點位數據。")

    except Exception as err:
        raise HTTPException(status_code=400, detail=f"檔案解析失敗: {str(err)}")

    sync_status = "本地解析完成（未填寫 Google Apps Script 網址，因此未同步至雲端）"

    if gas_url.strip() and gas_url.startswith("https://script.google.com/"):
        try:
            rows_payload = []
            for r in records:
                rows_payload.append([
                    patient["name"], patient["dob"], r["date_time"], r["zh_title"],
                    r["total_points"], r["avg_entropy"], r["entropy"]["level_1"],
                    r["entropy"]["level_2"], r["entropy"]["level_3"], r["entropy"]["level_4"],
                    r["entropy"]["level_5"], r["entropy"]["level_6"]
                ])
                
            payload = {"records": rows_payload}
            response = requests.post(gas_url.strip(), json=payload, timeout=20)
            
            if response.status_code == 200:
                resp_json = response.json()
                if resp_json.get("status") == "success":
                    sync_status = f"成功同步！已寫入 {len(rows_payload)} 筆檢測資料至試算表。"
                else:
                    sync_status = f"GAS 執行失敗: {resp_json.get('message')}"
            else:
                sync_status = f"GAS 請求異常，HTTP 狀態碼: {response.status_code}"
                
        except Exception as e:
            sync_status = f"呼叫 Google Apps Script 時發生錯誤: {str(e)}"

    return {
        "status": "success",
        "data": {
            "patient": patient,
            "records": records
        },
        "message": sync_status
    }
