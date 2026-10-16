# RGFCA 1,000種：2026-09-27版の真正な原本写真IDへの移行（2026-10-16）

## 2026-10-16 本番障害の確定と再取得しない復旧経路

[Source run #38847189688](https://github.com/zuizui0223/fcp/actions/runs/38847189688) は **failure** で終了（source diagnostic artifact [12029781544](https://github.com/zuizui0223/fcp/actions/runs/38847189688/artifacts/12029781544)）。元36GB dated archiveと4表のSHA256は成功したが、1,000対象種のうち**357種**で独立監査済み元写真数に対し抽出結果が不足し、全3,430,502 ID源CSVも128 workerも**未完成**。画像ピクセルの本番ROI測定数は **0件**。

実行runは `e2a23434ec171e3f2b7474c4972ab5cd91e31d36` の古い固定抽出コードをcheckoutしており、源写真IDが元42,111種対象の有効観察に紐づかない場合も保守的に衝突と判定した。これは元の全対象カウントの写真ID採用条件と一致しないため、357種の写真ID過剰隔離につながった可能性が高い。**推定原因であり、修正版の本物の3,430,502個別ID一致はまだ未実行。** PR現行ソースでは、全42,111候補のresearch-grade＋有効位置＋有効観察者を満たす外部観察だけが競合判定に寄与する。合成fixtureで完全原本importerとのPhotoID集合一致を確認済みだが、実データの成功と同一視しない。

原本36GBの再ダウンロードを繰り返さないため、main に[既存4表を再利用する手動専用Actions](https://github.com/zuizui0223/fcp/blob/main/.github/workflows/fcp-rgfca-1000-existing-monthly-tables-manual.yml) を追加した（GitHub commit `fee98a17168c6dc592d5d3519a64bf1ea7b256ef`）。これは `workflow_dispatch` のみで、研究者が事前に保管済み原本4表を準備した**専用永続ストレージ**と`fcp-rgfca-durable`ラベル付きself-hosted runnerでしか動かない。アーカイブを再取得するコマンドはなく、GitHub-hostedの失敗runをrerunするものでもない。

運用者の準備と手動順序：

- Repository Variables `FCP_RGFCA_ROOT`：一時ディスクではない永続パス。そこに `frozen_4_original_metadata_tables/{taxa,observations,observers,photos}.csv.gz` と `fcp_source_archive_pin.json`、凍結済み元42111種frame、独立監査count CSV+JSON、先験固定1000種roster+receiptを旧runの所定相対パスで配置する。これらは実在する検証済み原本ファイルが必要で、リンクや同名代用品は禁止。
- `FCP_RGFCA_EXISTING_PROVIDER_TABLES_APPROVED` に `I confirm these are existing approved persistent 2026-09-27 source tables; no archive redownload.` と記録し、実際の永続性・データ利用許諾も別途確認。**4表と別に200GBを要求しない**が、既存表を格納したうえで20GB以上の追加空きを技術上の事前条件とする。
- 手動 `preflight` で元frame/count/roster、独立audit SHA、4元表SHAと元archive pinを照合。次に `source_ids` で修正版コードが元3,430,502 IDを**全分類群ごとに完全照合**。結果が一致したら `screen_workers` で非結果依存の各64+64枚screenを選び128 blind workerのSHAを凍結。未完成/旧出力が既にある場合は上書きせず止める。

このmanual runの真正full-source成功artifactを、失敗済み `38847189688` のsuccessと偽って既存 `source→geo→pixel` のhardcoded run/SHAに渡さない。geoの新たな安全な起動接続は、sourceの実artifact・全count・128 SHAの検証が終わってから行う。**生物学的W/H2、広狭効果は現時点で未確定。**

## 今回の主作業と現在地

実データで確定済みの42,111種・93,209,351件は**植物写真の取得機会**であって花色測定ではない。事前に無作為抽出・固定した1,000種の対象原本機会数は **3,430,502件**。対象種名／種IDと200反復×250種の**種等重み**（各種50回）は変えていない。

**実行開始済み**：[FCP RGFCA EXACT 20260927 targeted 1000 real photo-ID source ONE SHOT — GitHub Actions #38847189688](https://github.com/zuizui0223/fcp/actions/runs/38847189688)。初回 #38846974428 は監査用の元ジョブが `main` で成功した事実を分岐の誤照合で拒否し、**原本取得前に失敗**。元SHAと成功状態を検証し直し、監査runだけ `main` を許可して修正した。これは過去の画像技術試験のlive APIではなく、SHA256を固定した**真正の2026年9月27日版 Open Dataset原本**を対象とする。ユーザーによる本番画像解析・多型種数への昇格はしない。

### 何を変更したか

- `scripts/acquisition/extract_rgfca_1000_original_monthly_photo_ids.py` — 元の4表（taxa、observations、observers、photos）のSHAを検査し、**全42,111種を対象とした種子孫の分類**を保持したまま、実際にデータベースへ保存するのは今回抽選した1,000種の観察・写真だけに変更。重複除去で生物学的判定を歪めないため、original photos表を**2回**読み、抽選対象外の観察で使われた写真IDとの重複は元42,111種の有効な観察同士に限って隔離する（無効な外部観察の再利用だけでは削除しない）。全種対象の独立監査済み写真数と**1,000種すべてで一致**しなければ、写真ID確定と宣言しない。撮影者と観察者を分離し、著作者名・ライセンスも保存する。途中経過・エラーは残す。
- `scripts/acquisition/plan_rgfca_screen_from_targeted_original_monthly_photo_ids.py` — 合格済み3,430,502件の元写真IDだけを入力に、**各種最大64枚×独立撮影者2組**のハッシュ・地理セル分散抽選を実施。42,111種の機会数台帳は残すが、未抽選種を「単色」としない。
- 既存の `scripts/acquisition/freeze_fcp_photo_first_blind_workers.py` への接続 — 新規原画ROI-v4/EfficientSAMへ公開するのは、`measurement_id,image_filename,photo_license` の3列だけ。種・座標・過去の色値を画像評価器へ渡さない。128分割を維持する。

### 検証の境界

[合成原本・完全経路の回帰検証（204 passed）](https://github.com/zuizui0223/fcp/actions/runs/38846581504) 成功。人工的に異なる種へ転用された写真IDは排除し、**原本の元種別写真数に一致しない場合は完全拒否**。観察のownerと写真のcreatorの違い、実写真の拡張子JPG、元著作者クレジット、原本種IDの変更拒否、結果への花色・文献の混入拒否を確認した。

成功すると one-shot ワークフローは `fcp-rgfca-1000-EXACT-20260927-original-photo-ID-128-worker-source-NO-PIXELS` のartifactへ**全1,000種の3,430,502写真ID台帳と128分割の予定写真ID一覧**を出す。ただし完了をチェックしない限りファイルは存在するとみなさない。失敗時は、相違する種の元写真数、どの段階で処理が途切れたかを `fcp-rgfca-1000-dated-source-INCOMPLETE-or-COLLISION-mismatch-NOT-ROI` の失敗証拠として保存し、**別のデータに置換して強制的に成功判定を出さない**。

## 次の生物学的条件

1. 真正な原本写真ID×元画像の一致と、本物の花領域のROI-v4測定完了。
2. 一次測定では撮影者・観察ID・空間セルと、花序内の白＋黄のような**構造由来色混在**を独立評価。
3. 1,000種×200反復の多型検出は `P_eval`、`P_poly|eval`、`P_joint`、**写真IDの実共有率**を別々に報告。96/128枚の反復では実際の写真集合が入れ替わることが必要。
4. 実写真12種（1,034枚）の過去の技術試験から認識したROI分類可能率（405/1,034=39.17%）は**本研究の普遍的な植物種間の確率とみなさない**。旧12種photo IDはlive API技術試験由来であり今回の正式母集団へ流用禁止。

**注意：この段階は原写真IDのsource台帳の取得であり、花色・W/H2・気候・土壌・遺伝的な多型の新たな結果ではない。** GitHub-hosted runnerでは時間／ディスク上限を超える可能性がある。成功・失敗の実行ログで実際に確定した範囲のみを採用する。
