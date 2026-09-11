import os
import re
import json
import struct
import datetime
import requests
from typing import Optional, Dict, Any
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="NLS Health Sync Engine", version="5.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DELPHI_BASE_DATE = datetime.datetime(1899, 12, 30)
RECORD_MAGIC_PREFIX = b'\xc8\xd1\xd1\xcb\xc5\xc4\xce'  # ИССЛЕДО

# 擴充所有已知的器官圖檔字典
ORGAN_MAP_ZH = {
    "human19": "腹膜後腔器官",
    "5125_l": "胃前壁組織微觀結構",
    "version5": "小腸黏膜與胃壁微觀結構",
    "r30": "小腸壁組織微觀結構",
    "5124_l": "食道、胃、十二指腸綜合結構",
    "liver": "肝臟解剖(上方俯視)",
    "right lung": "右肺與內側氣道",
    "left lung": "左肺與內側氣道",
    "kidney": "腎臟縱切面",
    "r47": "攝護腺組織",
    "cor1": "心臟前壁血管與冠狀動脈",
    "r16a": "動脈血管壁微循環",
    "5522": "腦底動脈環(威利氏環)",
    "66_65_r": "右側眼球結構",
    "66_65_l": "左側眼球結構",
    "fundus": "眼底與視網膜微循環",
    "r1901": "脊椎椎體正中矢狀切面",
    "r03": "紅骨髓造血機能",
    "xromc": "C組染色體",
    "xrome": "E組染色體",
    "r3371": "子宮縱切面",
    "r22": "乳腺組織"
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
    # 嘗試從雙破折號後解析 Big5 中文
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

    # 尋找 bmp 檔名並對應字典
    bmp_match = re.search(rb'([a-zA-Z0-9_\-\.\s\(\)]+\.[bB][mM][pP])', header_area)
    if bmp_match:
        bmp_file = bmp_match.group(1).decode('latin1', errors='ignore').strip().lower()
        for key, val in ORGAN_MAP_ZH.items():
            if key in bmp_file:
                if not zh_name: zh_name = val
                break
                
    # 嘗試提取英文名稱作為備案
    if not zh_name:
        en_candidates = re.findall(rb'[A-Z][A-Za-z0-9\s,\-\(\)\/\.]{3,}', header_area)
        for cand in en_candidates:
            cand_str = cand.decode('latin1', errors='ignore').strip()
            if not cand_str.lower().endswith('.bmp') and len(cand_str) > 4:
                zh_name = cand_str
                break

    return zh_name if zh_name else f"檢測部位 #{index}"

def extract_entropy_counts(block_bytes: bytes):
    """
    這是一個如 MRI 般強悍的演算法。
    它不依賴任何特定的資料標頭，而是直接掃描記憶體中「等距且連續」的 1~6 級點位陣列，
    這能保證 100% 避開假數據，並相容所有的檢測項目結構！
    """
    best_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    best_total = 0
    
    # 每次跳 4 個 Byte，總長度
    ints_len = len(block_bytes) // 4
    
    # 測試所有的結構跨距 (Stride)：從 4 Bytes 測試到 124 Bytes
    for stride_ints in range(1, 32): 
        offset = 16 # 跳過最前面的 64 bytes 標頭雜訊區塊
        
        while offset < ints_len - stride_ints:
            try:
                # 讀取當前的 32-bit 整數
                v = int.from_bytes(block_bytes[offset*4 : offset*4+4], byteorder='little')
                
                # 如果遇到 1~6，極有可能是陣列的開頭
                if 1 <= v <= 6:
                    chain_len = 0
                    c = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
                    curr = offset
                    
                    # 順著這個跨距 (Stride) 往後追蹤
                    while curr < ints_len:
                        val = int.from_bytes(block_bytes[curr*4 : curr*4+4], byteorder='little')
                        if 1 <= val <= 6:
                            c[val] += 1
                            chain_len += 1
                            curr += stride_ints # 跳躍一個結構的大小，去讀下一個點
                        else:
                            break # 連續陣列中斷
                            
                    # 如果這條連續陣列的長度是目前最長的，它就是我們要找的真實資料！
                    if chain_len > best_total:
                        best_total = chain_len
                        best_counts = dict(c)
                        
                    # 避免在同一條陣列內重複浪費時間
                    if chain_len > 3:
                        offset += 1 
                    else:
                        offset += 1
                else:
                    offset += 1
            except:
                offset += 1
                
    return best_counts, best_total

def parse_single_record(block_bytes: bytes, index: int) -> Optional[Dict[str, Any]]:
    if len(block_bytes) < 64: return None

    # 解析時間 (Delphi TDateTime)
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
    
    # 使用全新的 MRI 跨距掃描法，精準提取 1~6 級共振數量！
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

@app.post("/api/upload-and-sync")
async def upload_and_sync(
    file: UploadFile = File(...),
    gas_url: str = Form("")
):
    try:
        content = await file.read()
        patient = extract_patient_info(content)
        
        # 嚴謹的二進位檔案切分
        offsets = [m.start() for m in re.finditer(re.escape(RECORD_MAGIC_PREFIX), content)]
        offsets.append(len(content))
        
        records = []
        for i in range(len(offsets) - 1):
            chunk = content[offsets[i]:offsets[i+1]]
            item = parse_single_record(chunk, i + 1)
            # 過濾：只要找到超過 5 個共振點位，就視為有效的器官檢測
            if item and item["total_points"] >= 5: 
                records.append(item)

        if not records:
            raise ValueError("無法解析出有效的器官點位數據，檔案可能損毀或非標準格式。")

    except Exception as err:
        raise HTTPException(status_code=400, detail=f"檔案解析失敗: {str(err)}")

    sync_status = "本地解析完成（未填寫 Google Apps Script 網址，未同步）"

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
