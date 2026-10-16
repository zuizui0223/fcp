# FCP PR200 — 42,111種・実データ全写真ID凍結への引き継ぎ（2026-10-16）

## 新規：再取得なしの永続ストレージ手動Actions（2026-10-16）

mainに [`.github/workflows/fcp-pr200-durable-source-manual.yml`](https://github.com/zuizui0223/fcp/blob/main/.github/workflows/fcp-pr200-durable-source-manual.yml) を追加した。**workflow_dispatchのみ**、`push` / `pull_request` / `workflow_run` による自動起動はない。旧fallback失敗runを再実行せず、source側の修正済みPython module entrypoint（GitHub commit `feed51872e10611a7798f3db6c8c86a3b5d53527`）を固定checkoutする。GitHub-hosted runnerは使わず、専用ラベル `[self-hosted, linux, x64, fcp-durable]` を要求する。

**事前に運用者が準備するもの：**

- リポジトリ変数 `FCP_DURABLE_ROOT` に **/tmp でも /home/runner/work でもない、実際に永続なマウント配下の絶対パス**を登録する。専用Linux self-hosted runnerに `fcp-durable` ラベルを付ける。
- `$FCP_DURABLE_ROOT/FCP_DURABLE_STORAGE_APPROVED` に正確に `I confirm this is persistent FCP source storage, not hosted-runner scratch.` を記録する。**このファイルだけで物理永続性が証明されるわけではなく、運用者の確認を別途要する。**
- `inputs/42111_species.csv`（SHA256 `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`）、`references/all_42111_species_photo_opportunity_counts.csv`（`ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9`）、`references/all_42111_photo_opportunity_audit.json`（`ed62aee0f9730ad17347d25180575af1cd101e0fc35c581e470c4a17866ebf56`）を、**それぞれ成功済み原本count/独立audit成果物から**所定パスに置く。COUNT CSV/監査JSONの名前を誤って別ファイルで代用するとチェックで拒否される。
- `snapshots/inaturalist-open-data-20260927.tar.gz` は、**既に正当な取得権限の下で存在する固定日付原本**を使う。35,966,549,880バイト・SHA256 `bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e` が不一致なら停止。手動Actionsにはproviderからのアーカイブ取得処理は存在しない。原本を調達する許諾や永続保管がなければ**処理を開始しない**。
- `preflight` と `stage` は当該永続領域の**空き容量200,000,000,000バイト以上**、`emit` / `finalize` は60,000,000,000バイト以上が必要（総容量ではなく実空き容量）。依存Pythonパッケージはself-hostedに事前構築し、Actions中には無断の`pip install`も行わない。

**手動実行順序：** Actionsの `FCP PR200 durable 2026-09-27 original photo-ID source (manual)` で、まず `preflight`（原本/3参照SHA・200GB空き・運用者sentinel検査）、次に一度だけ `stage`（4表を走査して約9,320万の実IDを唯一のSQLiteとして永続保存）、その後 `emit` を `shard_start=0, shard_stop=4`、続いて `4,8`、`8,12`、… `124,128` と最大4区画ずつ計32回、最後に `finalize` を指定する。既存の完全seal区画はSHAを再検証し上書きしない。途中失敗で残った未seal区画は**手動原因検査と隔離が必要**で、自動で欠測を補完しない。

各runのartifactは**小さなsource stage/status/監査JSONとprogress traceだけ**。巨大なSQLite、36GB原本、全写真ID manifests、写真pixelsをGitHub Actions artifactとして転送しない。これらは永続領域にとどまる。

成功判定は `finalize` が `REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS` の独立照合receiptを出し、42,111全checkpoint・128 worker・93,209,351元photoIDと全分類群count一致を確認した時に限る。**この時点でROI測定済み0件であり、花色生物学的全件結果、W/H2、環境解析は未完了。** 初回本番64枚は、別途正式なsource admission、元ROI-v4/EfficientSAM fingerprint、提供者の取得許諾と全体レート管理を確認してから扱う。

## 実際に成功した工程（花色測定ではない）

- 原本：iNaturalist Open Dataset `metadata/inaturalist-open-data-20260927.tar.gz`（2026-09-27）。
- 原本 archive SHA256：`bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e`。
- 原本42,111種 frame SHA256：`5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`。
- [実計数 Actions #38749681793](https://github.com/zuizui0223/fcp/actions/runs/38749681793) は `success`、artifact `11992822641`。
- [独立監査 Actions #38818718877](https://github.com/zuizui0223/fcp/actions/runs/38818718877) は `success`、artifact `12012131967`。42,111件の元の種ID・学名を original frame と完全照合。
- [128分割の負荷計測 Actions #38819205800](https://github.com/zuizui0223/fcp/actions/runs/38819205800) は `success`、artifact `12012636357`。
- source photo IDs **93,209,351**件。写真1件以上 **42,095種**、0件 **16種**。40件以上 **32,794種**、100件以上 **26,510種**。
- 元のfcp固定 `FCP_ALL_OBSERVABLE_20261014` 種ハッシュによる128分割は変更していない。64枚/batchで **1,456,459 compute units**（写真・種のサンプリング上限ではない）。最大partition 005 = 1,260,277件、最小partition 025 = 405,702件。
- 厳密な原本count CSV SHA256：`ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9`。
- 本リポでの短い canonical provenance：`results/fcp_real_42111_source_opportunity_20261016/verified_source_counts_capacity_receipt.json`。

**生態学的な限界：** ここでいうphoto IDsは「原本の提供元から利用可能と判定された植物写真」であり、花を写した／ROI-v4で評価可能／白か有色かとは無関係。先行96実写真pilotは別日の実験集団であり、今回の9,320万件の画像測定数には決して足さない。種内多型W、H2、999-world null、土壌・気候の生態効果は**未測定**である。

## 次の唯一の本線：9,320万件の実写真ID universe を完全に凍結する

**推奨する実行基盤：** 十分な永続ディスク容量と実効的なジョブ再開機能を持つ機械。以前のGitHub hosted runnerでの全ID展開はディスク不足 `Errno 28`、続く復旧は約355分の実行上限で中断した。GitHub Actions artifactに小さなsource-count CSVしか残しておらず、元の36GBアーカイブと巨大なSQLite/写真ID manifestは**永続化されていない**。成功した「counts」を、本番の全写真ID凍結の代用にしないこと。

以下の`$FRAME`, `$ARCHIVE`, `$ROOT`, `$COUNT_CSV`, `$AUDIT_JSON`は永続ディスク内の真正な原本ファイルを指す。URIやdigestが違えば停止。既に取得した原本がある場合は再取得しない。適切な提供元バルクアクセス・利用規約、事前の保存/計算費用を確認した実行に限る。

```bash
# 1. 既に検証した同じ immutable dated snapshot と original species frame を使用。
# --stage-only は全写真IDをグローバル一意化したSQLiteに固定するが、
# この時点では42,111種すべての manifest / 画像pixels はまだ未完成。
python scripts/acquisition/execute_fcp_42111_census.py \
  --frame "$FRAME" --archive "$ARCHIVE" --work "$ROOT" \
  --expected-species 42111 --shards 128 --worker-partitions 128 \
  --stage-only --release-verified-temporary-tables

# 2. 同じ original source SQLite を再ダウンロードせず少数shardずつ生成。
# 各ブロックのあとにsource shard receiptを検証し、永続保存する。
python scripts/acquisition/emit_fcp_staged_source_shards.py \
  --frame "$FRAME" --staged-metadata "$ROOT/complete_42111_metadata_source" \
  --expected-species 42111 --shards 128 --shard-start 0 --shard-stop 4
# 以降 4:8, 8:12, ... 124:128 を独立した実行で継続し、未完了を成功扱いしない。

# 3. 全128shardを確定して、真の3列location-blind ROI workerを構築。
python scripts/acquisition/finalize_fcp_staged_source.py \
  --frame "$FRAME" --work "$ROOT" \
  --expected-species 42111 --shards 128 --worker-partitions 128

# 4. 今回の独立監査済み93,209,351件との全種・全count一致を確認。
# 成功しない限り原画像測定の本番母集団へ昇格させない。
python scripts/acquisition/verify_fcp_complete_photo_ids_vs_original_counts.py \
  --baseline-counts "$COUNT_CSV" --independent-audit "$AUDIT_JSON" \
  --complete-source-counts "$ROOT/all_42111_species_source_photo_counts.csv" \
  --full-source-receipt "$ROOT/full_source_execution_receipt.json" \
  --firewall-freeze "$ROOT/frozen_location_blind_all_photo_firewall/fresh_study_freeze.json" \
  --output "$ROOT/independent_original_all_ID_source_admission.json"
```

上記の最終監査で `REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS` が成立することが、**実写真のROI測定に移る条件**。

## 最終段階（原画像）

`scripts/acquisition/execute_fcp_full_photo_work_queue.py` を、既存の `docs/FCP_DISTRIBUTED_ALL_PHOTO_HANDOFF_20261015.md` どおり、元のROI-v4 YOLO/EfficientSAM SHAを確認してから実行する。64枚の単位は計算・リトライ単位であり、93,209,351件からの標本抽出や100枚／種制限を意味しない。

**提供元負荷の制約：** 現行codeのHTTPアクセス間隔は各workerで最低1.1秒で、別ホストの128worker全体にまたがる使用量は自動的には制御されない。独立した許諾・全体レート調整がなければ128個を並列公開URLダウンロードしない。全ID manifest・写真本番計測には永続保存/十分な演算資源が必要。短時間の小さなCIジョブで93,209,351件の花色分類済みと報告しない。

外部計算機・永続データ領域の準備はこのチャットのGitHub接続だけでは実行できないため、コード・原本SHA・検証法・128分割負荷表を先に確定し、未完了部分は明示する。
