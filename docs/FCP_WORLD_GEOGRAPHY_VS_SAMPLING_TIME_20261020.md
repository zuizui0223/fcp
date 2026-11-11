# FCP：世界花色の地域分化と観察年月の交絡（実測447,500枚）

**証拠更新 2026-10-20。** これは地理的写真色ラベル差を、観察年月の異なる地域構成と区別するための探索的解析であり、集団史・選択・色素分子機構の因果推定ではない。

## 元データと再現性

- [実測データの全工程成功 #39425739017](https://github.com/zuizui0223/fcp/actions/runs/39425739017)
- [完全な数値JSON](../results/fcp_world_geographic_photo_month_null_actual_20261020.json)
- [実装コード](../scripts/analysis/fcp_world_phenology_conditioned_geography_20261020.py)
- 元CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6
- 実測447,500枚のうち4色分類可能227,107枚、観察年月が2000–2026で有効226,601枚。日時不明・対象外は506枚。
- 通常画像と高純度画像の両方で、同一1,281種・同一2,762種×等面積格子。通常67,302写真／高純度55,806写真。
- 各種・格子の入場条件は10枚以上かつ3観察者以上。高純度は花のROIに割当色が65%以上・競合グループは20%以下。**高純度は遺伝的morph同定ではない。**

## 競合する3つの帰無モデル

各種の写真ラベルの総数、格子、観察者識別、写真・観察ID、観察年月と種の構成を保存する。全変動距離TVで同一種の格子間4色構成差を計算し、各種等重みで平均。

1. **unrestricted:** 同種内の任意の写真の色ラベルを交換。
2. **calendar_month:** 同じ暦月内だけ色を交換。
3. **calendar_year_month:** 同じ年・月内だけ色を交換。

各条件につき199回の色ラベル置換。主量は観測格子間TV−対応帰無期待。白色写真頻度差と非白色残差も計算するが、白色固有の適応優位とは解釈しない。

種内の同月・別地域に、観察者と観察IDが異なる写真ペアを5件以上持つ機会は、高純度で1,265種。同一年月・別地域で同じ条件なら1,083種であった。月–年月内の交換が可能な種を明示し、識別不可能な帰無を生態学的陰性結果と呼ばない。

## 主要実測結果

| 4色地域間TVの観測−null平均 | 通常写真（1,281種） | 高純度写真（同じ1,281種） |
|---|---:|---:|
| unrestricted | +0.008381 | +0.007263 |
| 暦月ごと | +0.007200 | +0.006351 |
| 同じ年・月ごと | +0.004784 | **+0.003515** |
| 年月条件後の片側置換検定 | p=.005 | **p=.005** |

高純度の年月条件での超過は無制約の半分程度である。**約51.6%の低下は帰無モデルを変更した場合の推定量の減少**であって、季節性の因果的説明率や分散分割ではない。暦月や年月は花齢、開花フェノロジー、撮影努力の違いも吸収する。

### 種群別の年・月条件付き高純度結果

| Species-disjoint source cohort | 種数 | 地域色差の超過 | 探索的片側p |
|---|---:|---:|---:|
| Discovery | 127 | +0.007754 | .015 |
| Validation | 135 | +0.000340 | .415 |
| Prospective | 155 | +0.000007 | .480 |
| Wave4 | 147 | +0.000835 | .330 |
| New U100 | 717 | +0.004669 | .005 |
| **New U100除外** | **564** | **+0.002047** | **.075** |

**判定：全プールでは観察年月を揃えた地域色差が残るが、厳格な品質基準下の種群間一般則は未支持。** New U100の寄与が大きい。New U100を除いた同じ564種を高純度ではなく通常写真で解析すると、同年月条件後の超過は+0.004386、p=.005であり、色純度条件への感度がある。

元のProspective群の事前固定性は連続色空間W仮説に対するものであり、今回の年–月地域TVの未開封追試ではない。全て探索的、最低到達Monte Carlo p=.005。

## 世界パターン、歴史と環境、分子機構をどう分離するか

**形質の生起（H1）:** 別の事前固定158種では、連続色空間の白―有彩色方向W=.517246、構成保存帰無中央値=.457143、p=.001が支持済み。しかし写真の白頻度・今回のTV分解はWとは別の統計量。

**地理的な配置（H2）:** 元の三種群では花色差と地理距離の相関がBIO5条件付きの小さな環境関連より強い。今回の時間条件後の残差は、さらに写真色組成が地理的に偏る候補を示す。ただし距離は遺伝的移住・浮動・氷期の集団史を識別しない。**歴史優位は未証明。**

**環境（H3）:** [PR #248](https://github.com/zuizui0223/fcp/pull/248)の1,574種・86,114枚では、位置平年BIO5の白写真との探索的な正の種内関連（平均+0.03233 log odds/℃、属CI 0.01786–0.04728）がある。一方、種別の留保地域への環境予測改善は得られていない。今回の地理TV残差を気候適応の証拠とはしない。

**多型維持（H4）:** 50 kmでの写真色組成と、同一交配集団で遺伝的に異なる色型が長期共存することを区別する。既存の5 km・同年月条件は不支持。個体追跡と繁殖成功が未測定。

**分子（H5）:** Silene littoreaのF3h1・Myb1aと花弁限定PAL／植物体全体WALは独立公開文献に基づく機構仮説。今回の447,500枚の個体にRNAやHPLCやゲノムのデータは連結していない。白化色素合成コストと適応度は未同定。

なお[PR #249](https://github.com/zuizui0223/fcp/pull/249)のcomponent-specific基準で比較すると、白型の**比例的な地理分化優位は不支持**（厳格174種、p=.220）。生起しやすい方向と特別に地理分化しやすい方向は別問題である。

## 論文用の許容Results（英語）

Regionally sampled conspecific flower-colour labels retained a small excess of between-region differentiation after conditioning label reassignment on the same observation year and month. In 1,281 taxa with matched regional sampling and a strict floral-colour ROI rule, mean observed-minus-null interregional total-variation distance was 0.00351 (199-permutation one-sided p=.005), compared with 0.00726 under unrestricted exchangeability. This signal was heterogeneous across disjoint source-species cohorts and no longer clearly supported after excluding the largest cohort (564 taxa, p=.075). Such photographic phenotype sorting does not identify genetic population differentiation, environmental selection or colonisation history.

## 次の主検証

1. 栽培・外来、複色花、花齢を独立監査し、自生かつ花器官確定の個体・地域に限定して再解析。
2. 地理距離をそろえ環境だけが異なる地点、環境をそろえ距離が異なる地点を比較。実年月の気象と土壌・送粉者相を追加し、観察者・系統・空間ブロックで検証。
3. 中立SNP/個体群史のある植物群のみ、真の遺伝的IBD/IBE・創始者効果・過去分布連結性を別レイヤーに加える。写真の距離相関を集団史と呼ばない。
4. 色素合成コスト・利益には同一遺伝背景のPAL/WAL対照と繁殖・生存適応度を要求する。

**main・原著の凍結H1/H2/W/50km結果は変更していない。** 今回の新たな年月条件検証はPR #249の作業ブランチで実行・固定する。
