# CLAUDE.md

## Working rules

- Be concise. Avoid conversational filler.
- Make targeted edits rather than rewriting entire files.
- Stop and ask when requirements are ambiguous instead of guessing.

## Project context

免疫螢光（IF）影像定量 pipeline。核心程式在 `ifquant/`，說明在 `README.md`。

輸入是 Olympus E-M5 Mark II 拍的 4608×3456 8-bit RGB JPEG，三張一組拍同一視野：
B channel = DAPI 細胞核、G channel = MyD88、R channel = iNOS。

## Constraints

- 影像資料在 `EOC3/`、`HMC3 activated marker IF GMy88 R INOX/`、`test/`，
  共約 296 MB，已 gitignore。不要 commit。
- 主指標是 `cyto_mean`（細胞質環），不是 `cell_mean`。
- 任何強度比較都必須先做 L1 線性化 + 除曝光時間。**不要直接比 8-bit 值** ——
  相機是自動曝光，同一批綠色影像的曝光時間在 0.625s 與 0.769s 之間跳動。
- 已知硬限制：沒有 secondary-only control，所以加性 offset 是自由參數。
  效應量只能報區間，不要報單一數字。
- 只有 1 個 field/group 時不做統計檢定（pseudoreplication）。
