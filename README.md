# bio-segmentation — 免疫螢光影像的每細胞強度定量

把螢光顯微鏡拍的照片，量成「每一顆細胞一個數字」，然後比較 control 與 treatment。

本專案處理的是一個不太理想的情況：影像是**消費級相機在自動模式下拍的 JPEG**，
而不是科學相機的 raw TIFF。整套 pipeline 的存在理由，就是把這種影像救回可用的程度，
並且誠實標示救不回來的部分。

---

## 1. 輸入影像是什麼

| 項目 | 內容 |
|---|---|
| 相機 | Olympus E-M5 Mark II（消費級 Bayer 相機掛在顯微鏡上） |
| 解析度 | 4608 × 3456 |
| 格式 | 8-bit RGB JPEG（有損壓縮） |
| ISO | 影像 EXIF 都是 1600，但機身設定是 **ISO AUTO**（見下） |
| 曝光 | **自動**（`ExposureProgram = 2`，機身在 P 模式） |
| 白平衡 | **自動**（MakerNote `WhiteBalance2 = Auto`，與 EXIF `WhiteBalance = 0` 一致）|

### 機身實際設定（由 Super Control Panel 照片確認）

| 項目 | 設定 | 影響 |
|---|---|---|
| 曝光模式 | P（程式自動） | 快門在 0.625–0.769 s 間跳動，必須除曝光時間 |
| ISO | **AUTO** | 相機可自行改增益。本批 EXIF 剛好都落在 1600，但這是結果不是保證 |
| 白平衡 | **自動**（以檔案為準）| ❌ 通道增益逐張浮動，事後無法校正 |
| 色彩空間 | **sRGB** | 確認 L1 用反 sRGB EOTF 是對的轉換函數 |
| 畫質 | LN / NORM | 8-bit 有損 JPEG |
| 對焦 | MF | 焦平面手動固定 |
| 閃燈 | 關閉，補償 ±0.0 | 純落射螢光 |
| 銳利度 / 對比 | S±0 / C±0 | 未額外加銳化或對比 |
| Picture Mode | **Custom** | ⚠️ 不是 Natural/Flat，tone curve 未知 |

> **ISO AUTO 的處理**：感光度對增益是線性的，所以 L1 除完曝光時間後再除 ISO 增益
> （`core.py: ISO_REF`）。以 ISO 1600 為基準，因此本批數值不變；但若之後某張的 ISO
> 跑掉，這一層會自動吸收掉，不會混進訊號差異裡。跑 `run_pair` 時 L0 會先印出整批的
> 曝光/ISO/白平衡範圍，任何漂移都會在 log 最上面看到。

### MakerNote 實測（`PB154766.JPG`，exiftool 13.59）

```
Gradation            : Normal; User-Selected
Picture Mode         : i-Enhance; 2
Color Space          : sRGB
White Balance 2      : Auto
Shading Compensation : Off
Noise Filter         : Standard
```

**以檔案為準，不要以機身選單為準。** 機身面板顯示白平衡是固定「晴天」、Picture Mode
是 Custom，但檔案記的是 Auto 與 i-Enhance。MakerNote 記的是**這張影像實際套用**的設定，
選單顯示的是**現在**的狀態，兩者不一致時一律以檔案為準。

| 欄位 | 結果 | 影響 |
|---|---|---|
| Shading Compensation | Off | ✅ 機身沒做周邊減光補正，不會和 L2 flat-field 重複校正 |
| Color Space | sRGB | ✅ 確認 L1 的反 EOTF 選對了函數族 |
| Gradation | **Normal**（非 Auto）| ✅ 色調曲線固定，不隨每張直方圖變動 |
| Picture Mode | **i-Enhance** | ❌ **場景自適應**，見下 |
| White Balance 2 | **Auto** | ❌ 通道增益逐張浮動，無法事後校正 |
| Noise Filter | Standard | ⚠️ 機身有做空間降噪，改變了逐像素統計 |

> **i-Enhance 是目前最嚴重的問題。** Gradation 是 Normal 本來是好消息 —— 但 Picture Mode
> 的**基底**是 i-Enhance，Olympus 這個模式會分析每張影像的內容，套用**因場景而異**的
> 對比與飽和度增強。也就是說我原本擔心 Gradation Auto 會造成的「逐張不同的非線性轉換」，
> 換成從 Picture Mode 這條路徑發生了。
>
> 這對本專案的主比較是**系統性混淆**，不只是雜訊：control 與 treatment 視野的影像內容
> 本來就不同（這正是要量的東西），因此兩組會被套上不同的增強曲線，而增強量與訊號強度
> 相關。反 sRGB EOTF 無法還原它，因為那個轉換逐張不同且未被記錄。
>
> **白平衡 Auto 是同一類問題的第二個來源**：R/G/B 增益逐張由相機自行決定，G 通道的
> 絕對值因此不可跨影像直接比較。

**結論修正方向**：這兩項都無法事後校正，所以跨影像的絕對強度比較比原本認定的還要更弱。
`cyto_mean` 的組間差異必須視為**半定量的方向性指標**，效應量區間應再放寬。原本的
「沒有 secondary-only control 所以 offset 是自由參數」仍然成立，現在是三個限制疊加。

### 從檔案讀回拍攝設定 — `python -m ifquant.exif_report <資料夾>`

```
python -m ifquant.exif_report "HMC3 activated marker IF GMy88 R INOX"
```

機身選單顯示的是**現在**的設定，影像檔記的是**當時**的設定 —— 白平衡那次矛盾就是這麼來的，
所以會影響 L1–L3 的設定一律從檔案讀回，不靠選單截圖。

標準 EXIF 由 Pillow 讀。但真正決定 pipeline 假設成不成立的欄位（尤其 **Gradation**）
在 Olympus MakerNote 裡，Pillow 不解析，需要 exiftool：

```
brew install exiftool           # macOS
apt install libimage-exiftool-perl   # Ubuntu
```

沒裝 exiftool 時標準 EXIF 那半仍會印，MakerNote 那半標為無法取得，不會用猜的補。
報告會標出整批不一致的欄位，並在 `Gradation = Auto` 時明確警告。

其中 `ShadingCompensation` 要特別留意：若機身開了內建周邊減光補正，它會和 L2 的
retrospective flat-field **重複校正**同一件事。

每個視野拍三張，中間換螢光濾片。三張都是 RGB 檔，但各自只有一個 channel 有訊號：

| 檔名範例 | 有訊號的 channel | 內容 | 在 pipeline 裡的角色 |
|---|---|---|---|
| PB154765 | **B**（藍） | DAPI 細胞核 | **找細胞在哪裡**（segmentation） |
| PB154766 | **G**（綠） | MyD88 | **要量的訊號**（measurement） |
| PB154767 | **R**（紅） | iNOS | 本次未使用 |

### 通道抽取規則

每張圖只取一個 channel，程式裡是 `ifquant/core.py` 的 `CH` dict。

> **不能直接量 RGB 影像的灰階值。** ImageJ 的 "Mean gray value" 在 RGB 影像上算的是
> 三個通道的亮度混合，不是純綠色訊號。這是最常見的錯誤。

### 本次實際用到的檔案

資料放在 `HMC3 activated marker IF GMy88 R INOX/`（**不在 repo 裡**，296 MB 已 gitignore）：

```
control    = PB154765.JPG (取 B) + PB154766.JPG (取 G)
treatment  = PB154768.JPG (取 B) + PB154769.JPG (取 G)

flat-field 由該資料夾全部 12 張綠色影像估出：
  766, 769, 796, 799, 802, 805, 808, 811, 814, 817, 820, 823
```

檔名對應規則：連續三個編號 = 同一個視野的 blue / green / red。
要重跑的話，把整個資料夾放回專案根目錄即可，路徑寫在 `ifquant/run_pair.py` 的 `SRC`。

---

## 2. 為什麼不能直接用 ImageJ 的 `Mean`

這是整套 pipeline 的動機。四個問題都是從這批影像實測出來的，不是理論顧慮。

### (1) 相機在自動曝光 — 最嚴重的問題

12 張綠色影像的 EXIF 曝光時間：

```
0.6250s  ×5
0.7692s  ×7      →  相差 23%
```

而我們要量的效應本身大約只有 1.5 倍。**自動曝光的方向還是反的**：
樣本越亮 → 相機曝光越短 → JPEG 看起來一樣亮。
它會主動抵銷你要量的訊號。

所以每張圖都必須**除以自己的曝光時間**，才能還原成「每秒收到多少光」。

### (2) JPEG 的 gamma

8-bit 像素值不是光子數，中間隔著 sRGB 的 tone curve（gamma ≈ 2.2）。
「A 比 B 亮 10%」這種比值運算在 gamma 空間裡不成立，必須先做反 gamma 還原成線性。

同一組資料，處理方式不同，結論完全不同：

| 算法 | control vs treatment 比值 |
|---|---|
| 直接比 8-bit 值（= ImageJ `Mean`） | **1.04** |
| 反 gamma 線性化後 | 1.10 |
| 線性化 + 除曝光時間 | **1.35** |

### (3) 照明不均（vignetting）是「乘性」的，不是「加性」的

中心亮、角落暗。線性空間下，角落只有中心的 **56%**（衰減 44%）。

照明不均的數學模型是 `影像 = 訊號 × 照明(x,y)`，所以正確的修正是**除法**（flat-field division）。
ImageJ 常用的 rolling ball（`Process ▸ Subtract Background`）是**減法**，
那是給加性背景（相機底噪、散射光）用的模型，用在照明不均上是錯的模型。

### (4) 細胞太密 + 訊號飽和

- **Confluent**：細胞鋪滿整個視野，幾乎沒有 cell-free 區域可以拿來當背景。
  「在空白處畫幾個框取平均」這招在這批影像上做不到。
- **飽和**：含飽和像素（值 = 255）的細胞比例，control **13.5%**、treatment **33.3%**。
  treatment 有三分之一的細胞打頂了，這是自動曝光的後遺症。

---

## 3. 做法：L1 – L5 五層

### 核心原則

> **在藍色（DAPI）找細胞，在綠色（MyD88）量訊號。**

絕不能用綠色自己去 threshold 找細胞 —— 那樣訊號弱的細胞會被整個漏掉，
剩下的都是亮細胞，平均值自然偏高。這是 selection bias。

### L1 線性化 — `core.py: srgb_to_linear`, `load_linear`

**白話**：把 JPEG 的 8-bit 值還原成「光的量」，再除以曝光時間。

```
8-bit 值 → 反 sRGB gamma → 除以 EXIF 曝光時間 → 除以 ISO 增益 → 「每秒的光」
```

之後所有數字的單位都是 **linear radiance / 秒**，不再是 8-bit 值。

⚠️ 這一層是近似的：色彩空間確認是 sRGB，但 Picture Mode 是 **Custom**，Olympus 仍會
套自己的 tone curve，不完全等於標準 sRGB EOTF。
這是這條救援路線最主要的殘留誤差。

### L2 照明修正 — `core.py: estimate_flatfield`, `apply_flatfield`

**白話**：用 12 張圖疊起來取中位數，就能看出「哪裡是鏡頭比較暗」，然後用除法攤平。

原理：每張圖的細胞長在不同位置，但照明不均永遠在同一個位置。
把 12 張圖的每個像素取中位數，細胞就被平均掉了，剩下的就是照明形狀。

兩個關鍵參數：
- 每張圖先除以自己的中位數 → 避免「某張圖比較亮」污染照明估計
- 平滑 σ = **500 px** → 必須遠大於一顆細胞（約 350 px），否則細胞紋理會被當成照明除掉

> 這裡踩過一個坑：一開始 σ 設 60 px，估出來的照明場是斑駁的，
> 等於把細胞紋理也除掉了。改成 500 px 後才是乾淨的碗狀。
> **驗證方法**：把 flat-field 畫出來看，應該是平滑的中心亮、四角暗，看不到任何細胞形狀。

### L3 背景扣除 — `core.py: offset_from_darkest`

**白話**：減掉相機底噪 + 非專一性螢光，讓「零」是真的零。

L2 修的是乘性的照明，L3 修的是加性的背景，兩者不同，**都要做**。

理想做法是用一張 **secondary-only（只打二抗、不打一抗）** 的影像當基準。
目前沒有這張，所以退而求其次用「視野最暗的 1 percentile」估。

⚠️ 這在 confluent 影像上會**高估**背景（最暗的像素其實是細胞之間的縫，不是真空白），
所以程式會輸出一份 **offset sensitivity sweep**，見第 5 節。

### L4 分割 — `segment.py`

分兩步：

**(a) 找細胞核**（`segment_nuclei_classical`）
在藍色 channel 上：高斯平滑 → 去除緩慢變化的背景 → Otsu 二值化 →
距離轉換 + seeded watershed 切開黏在一起的核。

實測核直徑中位數約 **193–205 px**。

> 這裡沒有用單純的 Otsu：這批影像的 DAPI 背景很亮（約 160），核只有 215，對比很低，
> 而且核佔了 31% 的面積，直接 Otsu 會偏掉。所以先減掉緩變背景再 threshold。

**(b) 長出細胞範圍**（`segment_cells_seeded`）
以每個核為種子，沿著綠色影像的**邊緣梯度**往外長，鄰居互相擋住。
比 Voronoi 好，因為 Voronoi 完全不看影像內容、純粹幾何等分。

`max_expand = 80 px` 這個上限**非常重要**：

> 沒有這個上限時，一顆位在平滑區域的細胞會一直擴張，
> 實測有一顆「細胞」吃掉了 **304 萬像素（全視野的 19%）**，
> 旁邊的細胞則完全分不到細胞質。加上上限後，
> 細胞直徑中位數 251 px vs 核 193 px，合理。

另外程式會先做 **channel registration**：藍綠是兩次獨立曝光、中間換過濾片，
視野會位移。control 這組實測位移 **dx = −16.8 px**，修正後可用細胞數從 72 增加到 95。

### L5 量測 — `measure.py: measure_cells`

每顆細胞量三個區間：

| 區間 | 定義 | 用途 |
|---|---|---|
| `nuc` | 細胞核 | 核內訊號 |
| `cell` | 整個細胞範圍 | 全細胞 |
| **`cyto`** | cell 減掉 nuc（細胞質環） | **主指標** |

**為什麼主指標是 `cyto_mean`**，有兩個理由：

1. MyD88 是 cytoplasmic adaptor protein，主要在細胞質，這才是生物學上關心的區間。
2. 綠色影像裡細胞核看起來是**亮的**，而紅色 iNOS 影像裡細胞核是**暗的**（典型細胞質 marker 的樣子）。
   實測 `corr(DAPI, 綠色) = +0.79`，偏高。這暗示 DAPI 的發射光譜可能漏進綠色濾片。
   **量細胞質環可以完全繞過這個疑慮。**

每個區間都給 `_mean`（濃度概念，跟細胞大小無關）和 `_integrated`（總量，跟細胞大小有關）。

QC 排除規則：
- `flag_saturated` — 細胞內有任何飽和像素
- `flag_area_outlier` — 面積在 1/99 percentile 之外
- `flag_no_cyto` — 細胞質環少於 2000 px（絕對值門檻，不是比例。
  用比例會優先砍掉擁擠區域的細胞，那正是不能有偏差的族群）
- 核碰到畫面邊緣的細胞直接移除

---

## 4. 怎麼跑

```bash
pip install -r requirements.txt
```

```bash
python -m ifquant.run_pair --backend classical
```

輸出在 `results/`：

| 檔案 | 內容 |
|---|---|
| `per_cell_classical.csv` | 每顆細胞一列，含三個區間的數值與所有 QC flag |
| `summary_classical.csv` | 四個指標的 control / treatment / 比值 |
| `offset_sweep_classical.csv` | offset 敏感度掃描（見下節） |
| `qc_classical.png` | 九宮格 QC 圖：分割邊界、flat-field、分布、敏感度 |
| `log_classical.txt` | 完整執行紀錄 |

**每次跑完一定要看 `qc_classical.png`**，特別是：
- 左下的 flat-field 應該是平滑的碗狀，看不到細胞形狀
- 中間兩張的橘線（細胞）/ 藍線（核）邊界要合理，不能有一顆吃掉一大片

---

## 5. 目前結果

單位：linear radiance / 秒。control 95 顆細胞、treatment 62 顆通過 QC。

| 指標 | control | treatment | T/C |
|---|---|---|---|
| `nuc_mean` | 0.5153 | 0.7994 | 1.55 |
| `cell_mean` | 0.4085 | 0.6444 | 1.58 |
| **`cyto_mean`**（主指標） | **0.3266** | **0.5415** | **1.66** |
| `cyto_integrated` | 16,956 | 31,897 | 1.88 |

對照：直接比 8-bit `Mean` 只會得到 **1.04**。

### 三個限制（重要，不要略過）

**(1) offset 無法識別 → 效應量只能給區間**

沒有 secondary-only control，加性背景就是一個自由參數：

```
假設的 offset    T/C
0.0000          1.33     ← 假設完全沒有加性背景
0.1833          1.47
0.3186          1.66     ← 程式目前用的估計值
0.4583          2.15     ← 再高 control 就變負值，物理上不可能
```

**誠實的結論是 T/C 落在 1.33 – 2.15 之間，不是 1.66。**

收斂的唯一方法：補拍一張 **secondary-only（只打二抗）** 影像，同樣的光路與曝光設定。
這張現在補拍仍然有效，不需要跟原樣本同一批。**這是投報率最高的一步。**

**(2) 飽和不對稱，方向對 treatment 不利**

treatment 有 33.3% 的細胞含飽和像素，control 只有 13.5%。
程式把它們排除了，等於優先砍掉 treatment 最亮的細胞
→ **1.66 是偏保守的估計**，真值應該更高。

**(3) n = 1 field / group，不能做統計**

目前 CSV 裡的 cell-level SD / SEM **只能當 QC 看**。
在單一視野的數百顆細胞上跑 t-test 是 pseudoreplication，p 值會假性地極小。

統計的複製單位是 coverslip / 獨立實驗，不是細胞。
正確做法：跑完所有 field → 每個 field 得到一個數字 → 在 replicate 層級做檢定。

---

## 6. 檔案結構

```
ifquant/
  core.py       L1 線性化、L2 flat-field、L3 offset
  segment.py    L4 分割（細胞核 / 細胞範圍）+ channel registration
  measure.py    L5 每細胞量測、QC flag、offset 敏感度掃描
  run_pair.py   跑一組 control/treatment，輸出 CSV + QC 圖
results/        輸出（已含一次執行結果）
requirements.txt
CLAUDE.md       給 Claude Code 的專案規則
```

`--backend cellpose` 的路徑也已實作（`segment_nuclei_cellpose` / `segment_cells_cellpose`，
使用 Cellpose-SAM），但這台機器沒有 GPU，CPU 推論尚未跑出結果。
目前所有數字都來自 classical backend。

---

## 7. 下一步需要的東西

1. **12 個 field 的 control / treatment 分組**，以及有幾個獨立的 biological replicate
   → 有了才能跑完整批並做統計
2. **一張 secondary-only 影像** → 把 T/C 從區間收斂成單一數字
3. `EOC3/*/RoiSet.zip`（手動框的 ROI）→ 可以拿來量化自動分割 vs 人工標註的一致性（IoU / F1）

### 給未來拍攝的建議

如果還有機會重拍，這三件事會讓上面大半的救援工作變成不必要：

- 相機切到 **M 模式**：固定曝光、**ISO 從 AUTO 改成固定值**
- 白平衡從 **Auto 改成固定預設**（晴天或自訂色溫皆可，重點是固定）—— 否則 R/G/B
  增益逐張浮動
- Picture Mode 從 **i-Enhance 改成 Natural**（Gradation 維持 Normal）—— i-Enhance 是
  場景自適應的，是目前最大的定量障礙
- **Noise Filter 關閉** —— 空間降噪會改變逐像素統計
- 存 **RAW（.ORF）** 而不是 JPEG —— 直接得到線性的 12-bit 資料，L1 整層可以省掉
- 拍一張 **secondary-only** 與一張 **flat-field 參考**（均勻螢光片）
