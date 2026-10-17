# FCP・画像から花色多型を判定する研究設計（2026-10-16）

## 研究上の決定

**多型かどうかの決定を文献一覧から画像に切り替える。** 42,111種は単に写真測定の候補であり、既知の多型種42,111種ではない。これまで作成した`analysis/fcp-literature-confirmed-local-polymorphism-20261016`の「文献でC=1と判定された26種限定」案は**本研究の生物学的選抜には使用しない**。文献のC/Sラベルは別途外的妥当性確認にのみ利用可能であり、画像スクリーニングのデータ・実行順序へ混入させない。

**母集団・原本**：2026-09-27 iNaturalist Open Datasetの42,111種・93,209,351件の適格写真ID（[source count #38749681793](https://github.com/zuizui0223/fcp/actions/runs/38749681793)、[全種・学名独立監査 #38818718877](https://github.com/zuizui0223/fcp/actions/runs/38818718877)）。この9,320万件は**植物写真の測定機会**であって、花が写った写真や、色が分類できた写真の数ではない。

## 判定法：写真の形質と実在の個体群多型を混同しない

1. バイト単位で凍結した元のROI-v4 YOLO→EfficientSAMの花領域と9色パレットを使用し、位置情報・種名をモデルに見せず `white/yellow_orange/red_pink/blue_purple` の4状態に分類する。ROI失敗、色型曖昧、通信失敗は混ぜず、独立した**分類不能**とする。写真全体の背景色や葉の色を花色と誤認した画像は採用しない。
2. 同一種の40枚以上が花色分類可能で、少なくとも2状態、少数側が10%以上、双方が各3件以上の独立観察に由来し、各状態が2名以上の観察者で見られる場合を**画像多型候補**とする。
3. observer-IDで決定的に2分したdiscover/holdoutの両方にそれぞれ20件以上の評価可能写真を要し、**同じ2状態が双方の分割で各2件以上の独立観察から繰り返し確認される**場合だけ`photo_supported_two_colour_states_observer_disjoint`と判定する。複数写真を持つ同一観察や、同一観察内で相互矛盾する色分類は独立した型の証拠としない。
4. 40枚未満、1色のみ、2型だが持続性の確認不足はすべて「多型未検出／画像からは未確定」である。**生物学的に単色固定・遺伝的多型なし、とみなさない**。画像の複数色が同じ自然集団内で共存するかは、GPS精度・時間・観察単位を使った後続の地理解析で検証し、50km近接だけで交配集団の同一性を断定しない。
5. 4色分類とDは探索的表現型の色違いであり、遺伝型や対立遺伝子の多型証明ではない。`W`は画像から選抜した種で再計算するが、元の固定白色軸・連続量ゲート・999帰無世界は変更しない。選抜に使った写真をそのまま独立な確認証拠とみなさず、可能なら**さらに別の観察者・写真IDを持つ推定用データ**を事前固定する。

## 元データ93,209,351件から第一次画像スクリーニングへ

92百万件全画像を先に測定する代わりに、**色を一枚も開かず**に原本写真IDと観察者・座標から第一次スクリーニング予定を決める：

- `scripts/acquisition/plan_fcp_photo_first_screen.py`が既存の固定原本SQLiteから全93,209,351件のIDを検査する。
- 各種で観察者単位の独立2組へ分け、各組最大64枚、計最大128枚／種を地理セル（5度）に分散させたハッシュ選抜で固定する。これは一次スクリーンの**計算単位**であって、種の採用上限や「128枚で多型なし」という陰性判断ではない。
- 一次選抜は最大約538万枚となる。実際の必要画像数は写真の利用可能性、観察者数、地理重複によって減り得る。実行後に正確な母数・地域・観察者の偏りを監査する。
- 判定不能となった種は追加独立写真を段階的に測定する。花色が多型らしい種だけに色のよい写真を探す、という都合のよい追加採取は許さない。未分類の種も母集団と未判定台帳に保持する。
- このスクリーンは**画像測定が終わっていない**。ソースID選抜はあくまで測定予定表で、元のfull-census`FRESH_BLIND_METADATA_FROZEN`とは分離しなければならない。

実行コードの順序：

```bash
# 永続ディスクの真正な2026-09-27ステージでのみ実行。
# 36GB原本のSQLite取り込みが未完なら、まず既存の
# scripts/acquisition/run_fcp_durable_source_phases.py --phase stage を実行。
python scripts/acquisition/plan_fcp_photo_first_screen.py \
  --staged-source /persistent/fcp/source_20260927/complete_42111_metadata_source \
  --original-frame /persistent/fcp/inputs/42111_species.csv \
  --source-counts /persistent/fcp/references/all_42111_species_photo_opportunity_counts.csv \
  --photos-per-fold 64 --max-per-observer 8 \
  --outdir /persistent/fcp/photo_first_screen_selection
```

出力は `photo_screen_original_metadata_selected.csv`、`photo_screen_all_species_opportunities.csv`、`photo_screen_source_selection_receipt.json`。サンプルは元画像を取得しておらず、画像による多型候補はまだ0件も確定していない（0候補という生物学的結果ではなく、測定未実施）。

## 一次スクリーニングから元画像計測までの封印された実行経路

`photo_screen_original_metadata_selected.csv`ができたら、文献のC/S分類表には一切触れず、画像ワーカー入力を固定する。

```bash
# 1. 選抜された写真IDだけを128分割し、元の測定モデルには3列のみ公開。
python scripts/acquisition/freeze_fcp_photo_first_blind_workers.py \
  --selected-metadata /persistent/fcp/photo_first_screen_selection/photo_screen_original_metadata_selected.csv \
  --selection-receipt /persistent/fcp/photo_first_screen_selection/photo_screen_source_selection_receipt.json \
  --source-frame /persistent/fcp/inputs/42111_species.csv \
  --source-counts /persistent/fcp/references/all_42111_species_photo_opportunity_counts.csv \
  --outdir /persistent/fcp/photo_first_screen_firewall

# 2. 同じ元ROI-v4 YOLO・EfficientSAMモデルのSHAを検証してから、
#    選抜された64枚の1バッチだけを処理。実写真へのprovider-rate管理が必須。
export FCP_PROVIDER_GLOBAL_RATE_FILE=/persistent/fcp/provider_s3_rate.lock
bash scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --full \
  --sam-dir /persistent/fcp/original_sam
python scripts/acquisition/run_fcp_photo_first_roi_work_unit.py \
  --firewall /persistent/fcp/photo_first_screen_firewall \
  --partition 0 --batch 0 --batch-size 64 \
  --original-measure-script scripts/analysis/measure_random_photo_first_blind_partition.py \
  --original-detector data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt \
  --original-efficient-sam /persistent/fcp/original_sam \
  --output /persistent/fcp/photo_first_roi_terminals/partition_000

# 3. 全128分割について、存在する全バッチを重複なく処理した後だけ、
#    元photo ID↔blind測定IDを再結合する。
python scripts/acquisition/reassemble_fcp_photo_first_screen.py \
  --screen-firewall /persistent/fcp/photo_first_screen_firewall \
  --screen-terminal-root /persistent/fcp/photo_first_roi_terminals \
  --outdir /persistent/fcp/photo_first_roi_joined
```

ここで初めて`newly_remeasured_photo_screen.csv`と`sealed_photo_screen_terminal.json`が成立する。どのワーカーでも旧149,900枚の色結果を読んではならず、`run_fcp_photo_first_roi_work_unit.py`が使う画像評価器は凍結された元のROI-v4である。未完バッチ、偽のfull-census結果、異なるモデルを混ぜた結果は結合で拒否する。**このコマンド群はコードとして実装したが、永続実行環境上の原写真計測がまだ成功していないので、新規多型種数は未算出。**

## 分類済み写真だけから本解析対象を決める

`fcp_pipeline/photo_polymorphism_screen.py` は、完全な本研究ROI画像測定台帳または独立に封印した非網羅スクリーニングの計測結果だけを入力とし、元42,111種ID・元写真計数・元測定ハッシュと照合する。文献リスト、気候、土壌、遺伝情報から種を選ばない。

```bash
python -m fcp_pipeline.photo_polymorphism_screen \
  --measured /persistent/fcp/roi/photo_screen_new_measurements.csv \
  --measurement-receipt /persistent/fcp/roi/sealed_photo_screen_terminal.json \
  --source-frame /persistent/fcp/inputs/42111_species.csv \
  --source-counts /persistent/fcp/references/all_42111_species_photo_opportunity_counts.csv \
  --outdir /persistent/fcp/photo_first_polymorphism_cohort
```

結果は`photo_supported_polymorphic_species.csv`（画像判定が支持）、`not_proven_polymorphic_do_not_call_monomorphic.csv`（追加測定・品質審査対象）、`all_42111_photo_polymorphism_screen.csv`（全種台帳）。封印・結合のコードは実装済みだが、原写真の本番計測ジョブは未完了である。したがって`--measured`の生物学的実データと有効な`--measurement-receipt`はまだ得られておらず、**実写真から多型種が何種見つかったかは未算出**。

## 論文で説明できる範囲

本画像選択で問えるのは「画像から再現可能な複数花色型を認めた植物種において、その色差の量・白色軸・地域分割がどの環境異質性と共変するか」。検出率を世界の植物全体のFCP普及率とは呼ばない。非検出＝単色の対照群にしない。写真から遺伝的多型、色素の適応的利益、送粉者による自然選択を実証したとは主張しない。
