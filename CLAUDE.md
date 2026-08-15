# CLAUDE.md

## Working rules

- Be concise. Avoid conversational filler.
- Make targeted edits rather than rewriting entire files.
- Stop and ask when requirements are ambiguous instead of guessing.

## Project context

免疫螢光（IF）影像定量 pipeline。核心程式在 `ifquant/`，說明在 `README.md`。

輸入是 Olympus E-M5 Mark II 拍的 4608×3456 8-bit RGB JPEG，三張一組拍同一視野：
B channel = DAPI 細胞核、G channel = MyD88、R channel = iNOS。

## Tooling

`gh` CLI 已安裝並授權（帳號 anika826283，scopes: repo/workflow/read:org/gist），
但**不在這個 shell 的 PATH 上**。要用完整路徑呼叫：

```
$gh = "$env:LOCALAPPDATA\Microsoft\WinGet\Packages\GitHub.cli_Microsoft.Winget.Source_8wekyb3d8bbwe\bin\gh.exe"
```

PR / issue / CI 一律走 `gh`，**不要用瀏覽器自動化操作 GitHub** —— 先前試過，
`form_input` 設 PR 標題會被 React 狀態覆蓋，且編輯器的自動列表接續會弄亂 markdown。

## Documentation rule

任何分析只要產生結論，就要在 `docs/` 存一份 md，並在 `docs/README.md` 的索引表加一列。
如果後續開發或新資料**改變了既有結論**，回去更新對應的 md 與 `README.md`，
並在檔案開頭的「狀態」欄註明被推翻/被修正。不要只在對話裡講。

## Constraints

- 影像資料在 `EOC3/`、`HMC3 activated marker IF GMy88 R INOX/`、`test/`，
  共約 296 MB，已 gitignore。不要 commit。
- 主指標是 `cyto_mean`（細胞質環），不是 `cell_mean`。
- 任何強度比較都必須先做 L1 線性化 + 除曝光時間。**不要直接比 8-bit 值** ——
  相機是自動曝光，同一批綠色影像的曝光時間在 0.625s 與 0.769s 之間跳動。
- 已知硬限制：沒有 secondary-only control，所以加性 offset 是自由參數。
  效應量只能報區間，不要報單一數字。
- 只有 1 個 field/group 時不做統計檢定（pseudoreplication）。
- **control/treatment 分組目前未知。** 先前沿用的 `766=control / 769=treatment`
  已被時間戳推翻（相隔 21 秒）。不要在使用者確認分組之前報告任何效應量。
- **自動曝光混淆**：`cyto_mean = 線性化值 / 曝光時間`，而 AE 依整個視野的總亮度
  決定曝光，所以細胞密度會經由 AE 直接灌水 per-cell 數值。跨 field 比較時，
  一律要同時報告「扣掉曝光比值後的殘餘」，見 `ifquant/confound_check.py`。
