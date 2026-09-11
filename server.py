import os
import re
import json
import struct
import datetime
import requests
from typing import Optional, Dict, Any
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="NLS Health Sync Engine", version="4.0.0")

# 允許跨域請求，讓前端網頁可以呼叫 API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DELPHI_BASE_DATE = datetime.datetime(1899, 12, 30)
RECORD_MAGIC_PREFIX = b'\xc8\xd1\xd1\xcb\xc5\xc4\xce'  # 區塊開頭: ИССЛЕДОВАНИЕ

ORGAN_MAP_ZH = {
    "human19.bmp": "腹膜後腔器官",
    "5125_l.bmp": "胃前壁組織微觀",
    "version5.bmp": "小腸黏膜與胃壁",
    "r30.bmp": "小腸壁組織微觀",
    "5124_l.bmp": "食道、胃、十二指腸",
    "the liver, viewed from above.bmp": "肝臟解剖(上方)",
    "right lung (medial view).bmp": "右肺與內側氣道",
    "left lung (medial view).bmp": "左肺與內側氣道",
    "the kidney (left).bmp": "左腎縱切面",
    "the kidney (right).bmp": "右腎縱切面",
    "r47.bmp": "攝護腺組織",
    "cor1.bmp": "心臟前壁血管與冠狀動脈",
    "r16a.bmp": "動脈血管壁微循環",
    "5522.bmp": "腦底動脈環",
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
    """提取受檢者姓名與基本資料"""
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

def extract_binary_entropy(block_bytes: bytes) -> Dict[int, int]:
    """無差別滑動視窗：掃描最長連續 1~6 的整數陣列 (解決漏抓項目與假點位問題)"""
    best_overall_seq = []
    
    # 測試 4 種不同的位元組對齊 (Offset 0~3)
    for offset in range(4):
        current_seq = []
        best_seq = []
        # 以 4 bytes (32-bit Little Endian) 為單位進行底層掃描
        for i in range(offset, len(block_bytes) - 3, 4):
            val = int.from_bytes(block_bytes[i:i+4], byteorder='little')
            if 1 <= val <= 6:
                current_seq.append(val)
            else:
                if len(current_seq) > len(best_seq):
                    best_seq = current_seq
                current_seq = []
        if len(current_seq) > len(best_seq):
            best_seq = current_seq
            
        if len(best_seq) > len(best_overall_seq):
            best_overall_seq = best_seq

    counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    
    # 如果找到的連續陣列長度大於 5，這絕對就是真實點位資料，排除所有雜訊
    if len(best_overall_seq) >= 5:
        for v in best_overall_seq:
            counts[v] += 1
            
    return counts

def parse_single_record(block_bytes: bytes, index: int) -> Optional[Dict[str, Any]]:
    """解析單一器官區塊"""
    if len(block_bytes) < 64: return None

    # 1. 取得時間
    record_time = "未知"
    for offset in range(8, min(120, len(block_bytes) - 8)):
        try:
            val = struct.unpack('<d', block_bytes[offset:offset+8])[0]
            if 40000.0 <= val <= 55000.0: 
                dt = DELPHI_BASE_DATE + datetime.timedelta(days=val)
                record_time = dt.strftime("%Y-%m-%d %H:%M")
                break
        except: continue

    # 2. 取得名稱
    header_area = block_bytes[:1000]
    zh_name = ""
    
    bmp_match = re.search(rb'([a-zA-Z0-9_\-\.\s\(\)]+\.[bB][mM][pP])', header_area)
    if bmp_match:
        bmp_file = bmp_match.group(1).decode('latin1', errors='ignore').strip().lower()
        bmp_file = re.sub(r'^_+', '', bmp_file)
        for k, v in ORGAN_MAP_ZH.items():
            if k in bmp_file or bmp_file in k:
                zh_name = v
                break

    if not zh_name:
        dash_pos = header_area.find(b'--') 
        if dash_pos != -1:
            tail = header_area[dash_pos+2:dash_pos+160]
            tokens = []
            cur_pos = 0
            while cur_pos < len(tail) - 1:
                try:
                    ch = tail[cur_pos:cur_pos+2].decode('big5')
                    if '\u4e00' <= ch <= '\u9fff' or ch in '()（）-，、':
                        tokens.append(ch)
                        cur_pos += 2
                    else:
                        if tokens: break
                        cur_pos += 1
                except:
                    if tokens: break
                    cur_pos += 1
            if tokens: zh_name = "".join(tokens).strip()

    if not zh_name: 
        zh_name = f"檢測項目 #{index}"

    # 3. 取得 1~6 級點位
    counts = extract_binary_entropy(block_bytes)
    total_measured_points = sum(counts.values())
    
    weighted_sum = sum(k * v for k, v in counts.items())
    avg_entropy = round(weighted_sum / total_measured_points, 2) if total_measured_points > 0 else 2.5

    return {
        "index": index,
        "date_time": record_time,
        "zh_title": zh_name,
        "total_points": total_measured_points,
        "avg_entropy": avg_entropy,
        "entropy": {
            "level_1": counts[1], "level_2": counts[2], "level_3": counts[3],
            "level_4": counts[4], "level_5": counts[5], "level_6": counts[6]
        }
    }

@app.post("/api/upload-and-sync")
async def upload_and_sync(
    file: UploadFile = File(...),
    gas_url: str = Form("") 
):
    try:
        content = await file.read()
        patient = extract_patient_info(content)
        
        # 放棄正則運算，直接使用記憶體切割法，速度更快且絕不漏抓
        chunks = content.split(RECORD_MAGIC_PREFIX)
        
        records = []
        for i in range(1, len(chunks)):
            chunk = chunks[i]
            item = parse_single_record(chunk, len(records) + 1)
            # 過濾掉少於 5 個點的無效區塊
            if item and item["total_points"] > 5: 
                records.append(item)

        if not records:
            raise ValueError("無法解析出有效的器官點位數據，請確認檔案格式是否正確。")

    except Exception as err:
        raise HTTPException(status_code=400, detail=str(err))

    sync_status = "本地解析完成（未填寫 GAS 網址，因此未同步至雲端）"

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
            # allow_redirects=True 會自動跟隨 Google Apps Script 的 302 重新導向
            response = requests.post(gas_url.strip(), json=payload, timeout=20, allow_redirects=True)
            
            # 因為 GAS 會重新導向，只要狀態碼是 200 就是成功執行了 doPost (代表資料已寫入)
            if response.status_code == 200:
                sync_status = f"成功同步！已寫入 {len(rows_payload)} 筆檢測資料。"
            else:
                sync_status = f"同步可能有狀況，GAS 回傳狀態碼: {response.status_code}"
                
        except Exception as e:
            sync_status = f"資料已解析，但呼叫 GAS 時發生錯誤: {str(e)}"

    return {
        "status": "success",
        "data": {
            "patient": patient,
            "records": records
        },
        "message": sync_status
    }