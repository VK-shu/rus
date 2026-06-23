# RUS 失敗路徑明確化 — 實驗紀錄

> **日期：** 2026-06-23  
> **專案：** `SliQSim/external/rus`  
> **相關程式：** `rus/rus.cpp`、`es/exactdecomposer.cpp`

---

## 1. 實驗動機

修復 segfault 後，`rusSyn` 在部分參數下雖不再崩潰，但仍可能：

- PSLQ 找不到任何候選
- 某個 `OmegaRing` 在 Normalization 或 Decomposition 失敗
- 最終輸出**空 QASM** 或無意義結果，且**沒有清楚錯誤訊息**

為了方便實驗除錯與紀錄，將失敗路徑改為**明確回報**，不再靜默失敗。

---

## 2. 修改內容摘要

### 2.1 `exactDecomposer`

| 修改前 | 修改後 |
|--------|--------|
| `decompose(matr, c)` 為 `void`，中間找不到 generator 時 `return`，可能留下空電路 | 回傳 `bool`：`true` 成功、`false` 失敗 |
| `slC.find` 找不到查表結果時仍 `push_back` 空子電路 | 若 lookup 結果為空，回傳 `false` |
| 靜態 `decompose(matr)` 永遠回傳 `circuit` | 失敗時 `throw std::runtime_error("exactDecomposer: decomposition failed")` |

### 2.2 `RUS::run()`

新增以下檢查與行為：

| 階段 | 條件 | 行為 |
|------|------|------|
| Stage 1 結束 | `results.empty()` | `throw`: `RUS: PSLQ produced no candidates (...)` |
| Stage 2 | Normalization 拋例外（如 `Diverge`） | 記錄候選編號與原因，`continue` 下一候選 |
| Stage 3 | `decompose` 回傳 `false` | 記錄候選編號，`continue` 下一候選 |
| 全部候選處理完 | `success_count == 0` | `throw`: `RUS: all candidates failed normalization or decomposition` |
| 選優結束 | `best_cir.empty()` | `throw`: `RUS: no circuit selected after synthesis` |

`-D 1` 以上時，單一候選的 Stage 2/3 失敗會印到 **stderr**，程式會繼續嘗試其他候選。

---

## 3. 錯誤訊息對照

### 3.1 終止整次合成（exit code 1）

```
Error: RUS: PSLQ produced no candidates (try increasing -F, -I, or relaxing -E)
```

**常見原因：** `-I` 太小、`-F` 太小、`-E` 太嚴、角度過難。

---

```
Error: RUS: all candidates failed normalization or decomposition
```

**常見原因：** PSLQ 有候選但每個都在 Stage 2 或 Stage 3 失敗（參數邊界、數值問題）。

---

```
Error: RUS: no circuit selected after synthesis
```

**理論上不應出現**（防禦性檢查）；若出現代表選優邏輯異常。

---

```
Error: exactDecomposer: decomposition failed
```

**出現時機：** 其他模組直接呼叫靜態 `exactDecomposer::decompose()` 且分解失敗（非 `RUS::run` 路徑）。

---

### 3.2 單一候選失敗（不中斷，`-D ≥ 1` 可見）

```
RUS: normalization failed for candidate 3: Diverge
RUS: decomposition failed for candidate 5
```

---

## 4. 測試紀錄

### 4.1 正常路徑（應成功）

```bash
./rusSyn -T pi/16 -O out.qasm -F 5 -E 1e-5 -C g-count -D 1
```

**預期：** 寫出 `output/out.qasm`，exit 0。

---

### 4.2 PSLQ 無候選（應明確失敗）

```bash
./rusSyn -T pi/8 -O out.qasm -F 1 -E 1e-10 -I 1000 -D 1
```

**預期：**  
`Error: RUS: PSLQ produced no candidates (...)`  
不產生有效 QASM（或根本不寫檔，視 main 是否在 `run()` 後才寫檔）。

---

### 4.3 部分候選失敗、其餘成功（應成功）

```bash
./rusSyn -T pi/8 -O out.qasm -F 20 -E 1e-3 -C g-count -D 1
```

**預期：** 可能看到若干 `decomposition failed for candidate N`，但只要有一個候選成功即輸出 QASM。

---

## 5. 與先前 segfault 修復的關係

| 項目 | segfault 修復 | 本次失敗路徑 |
|------|---------------|--------------|
| 目的 | 避免 crash、確保找得到候選 | 失敗時給可讀錯誤、不產生空結果 |
| PSLQ 預設迭代 | 未給 `-I` → 1e8 | 不變 |
| `-I` 可選覆寫 | 支援 | 不變 |
| 多候選容錯 | 依賴多個 OmegaRing | 單候選失敗可跳過，全部失敗才報錯 |

詳見 [RUS_SEGFAULT_FIX.md](RUS_SEGFAULT_FIX.md)。

---

## 6. 後續可改進（未實作）

- 將失敗統計寫入 JSON/log 檔供批次實驗分析

### 6.1 已實作（2026-06-23）

- **選優邏輯**：`success_count == 0` 取代 `r_count == 1`，第一個**成功**候選才設為 `best_cir`
- **`count_t` / `count_gate`**：直接遍歷 `circuit` gate id，不再 serialize 字串計數
- **`normalization.cpp`**：`log2` / `ceil` 改為 `mpf_class` + MPFR，移除 `.get_d()` 路徑

---

## 7. 變更檔案

| 檔案 | 變更 |
|------|------|
| `es/exactdecomposer.h` | `decompose` 改回傳 `bool` |
| `es/exactdecomposer.cpp` | 失敗回傳 `false`；靜態版 throw |
| `rus/rus.cpp` | PSLQ / 逐候選 / 最終結果檢查 |
| `docs/RUS_FAILURE_HANDLING.md` | 本實驗紀錄 |

---

## 8. effort 預設與 CLI 輸入驗證（補充）

> **日期：** 2026-06-23  
> **相關程式：** `apps/main_rus.cpp`

### 8.1 問題

原先 `-F`（effort）預設為 `0`，PSLQ 外層迴圈 `while (_effort-- > 0)` **一次都不執行**，等於沒跑 PSLQ 卻不易察覺。  
其餘參數（`-E`、`-I`、`-C`、`-D`）也缺少啟動前檢查，錯誤值可能拖到演算法中途才失敗。

### 8.2 修改內容

**effort 預設改為 1**

```cpp
("effort,F", po::value<int>(&effort)->default_value(1), "Effort level (>=1). Default: 1")
```

未指定 `-F` 時，至少會嘗試收集 1 個 PSLQ 成功候選。

**新增 `validate_cli_params()`**（`main` 解析參數後、進入合成前呼叫）：

| 參數 | 條件 | 錯誤訊息 |
|------|------|----------|
| `-F` effort | `>= 1` | `effort (-F) must be >= 1` |
| `-E` epsilon | `> 0` | `epsilon (-E) must be > 0` |
| `-D` debuglevel | `>= 0` | `debuglevel (-D) must be >= 0` |
| `-I` iterations | 有指定時 `> 0` | `iterations (-I) must be > 0` |
| `-C` criterion | `g-count` 或 `t-count` | `criterion (-C) must be g-count or t-count` |
| `-P` precision | `--database` 時 `> 0` | `precision (-P) must be > 0` |

驗證失敗時由既有 `catch` 印出 `Error: ...` 並 **exit 1**，不進入 `RUS::run()`。

### 8.3 測試紀錄

**非法 effort：**

```bash
./rusSyn -T pi/16 -F 0 -E 1e-5
# Error: effort (-F) must be >= 1
```

**未指定 `-F`（新預設）：**

```bash
./rusSyn -T pi/16 -O out.qasm -E 1e-5 -C g-count -D 1
# 等同 -F 1，應正常產生 output/out.qasm
```

**非法 criterion：**

```bash
./rusSyn -T pi/16 -C foo -E 1e-5 -F 1
# Error: criterion (-C) must be g-count or t-count
```

### 8.4 變更檔案

| 檔案 | 變更 |
|------|------|
| `apps/main_rus.cpp` | effort 預設 1、`validate_cli_params()` |

---

## 9. PSLQ 精度、漸進迭代與 effort 預算（補充）

> **日期：** 2026-06-23  
> **相關程式：** `rus/rus.cpp`、`rus/rus.h`、`rus/pslq.cpp`、`apps/main_rus.cpp`

### 9.1 誤差判斷改為 mpf_class 高精度

**問題：** `gen_error_function` 與 debug 輸出使用 `toComplex(0)` 與 `double` 的 `cos`/`sin`，與建構子內 2048-bit MPFR 設定不一致。

**修改：**

- 新增 `omega_ring_to_mpf_complex(z, de)`：依 `toHprComplex` 同一公式，係數與 `cos`/`sin` 皆用 `mpf_class`
- 新增 `compute_phase_error(z, theta)`：PSLQ 成功判斷與 `-D 2` 的 `real epsilon` 皆走此函式
- 移除 `gen_error_function` 內的 `double` 路徑

### 9.2 PSLQ 漸進迭代（未指定 `-I` 時）

| 常數 | 值 | 用途 |
|------|-----|------|
| `INITIAL_PSLQ_MAX_ITER` | 1,000,000 | 第一次嘗試上限 |
| `DEFAULT_PSLQ_MAX_ITER` | 100,000,000 | 漸進加倍的上限 |

**流程：**

1. 以 1e6 次上限跑 PSLQ
2. 若 `r_sets` 為空，上限加倍（2e6 → 4e6 → …），直到有候選或達 1e8
3. `-D 1` 時印：`PSLQ: no candidates at limit N, retrying with M`

**指定 `-I N` 時：** 只跑一輪、上限為 N，不漸進加倍（與先前「CLI 控制迭代」行為一致）。

### 9.3 每輪 effort 重置 `it_count`

**修改前：** `it_count` 只在 PSLQ 開始時設為 `itm`，多輪 `-F` **共用**同一迭代預算。

**修改後：** 每收集一個成功候選後，除重置 `status` 外，另執行 `it_count = itm`，**每輪 effort 各自享有完整 `-I`（或漸進上限）次迭代**。

**語意：**

- `-F`：最多收集幾個成功候選（外層）
- `-I` / 漸進上限：每一個候選搜尋時，內層最多跑多少次 `pslq_iteration`

兩者不再共享同一個遞減計數器。

### 9.4 測試紀錄

```bash
# 預設漸進：多數角度 1e6 即足夠，不必每次等 1e8
./rusSyn -T pi/16 -O out.qasm -F 5 -E 1e-5 -C g-count -D 1

# 強制單次上限（不漸進）
./rusSyn -T pi/16 -O out.qasm -F 5 -E 1e-5 -I 500000 -D 1
```

### 9.5 變更檔案

| 檔案 | 變更 |
|------|------|
| `rus/rus.h` | `INITIAL_PSLQ_MAX_ITER`、phase error 宣告 |
| `rus/rus.cpp` | 高精度誤差、漸進 PSLQ 迴圈 |
| `rus/pslq.cpp` | 每輪 effort 重置 `it_count` |
| `apps/main_rus.cpp` | `-I` 說明文字更新 |
