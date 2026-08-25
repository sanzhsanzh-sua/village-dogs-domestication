# Village dogs / wild canids pilot — сводка находок

Пилотное сравнение village dogs vs wild canids (wolves, coyote) на 246
кандидатных доместикационных регионах (CDR) из Pendleton et al. 2018
(*Comparison of village dog and wolf genomes highlights the role of the
neural crest in dog domestication*, BMC Biology, PMC6022502).

## 1. Данные

Изначально запрошенный BioProject **PRJNA448733** оказался проектом
секвенирования **чистокровных пород**, а не village dogs (сверено с
Supplementary Data 1 статьи Plassais et al. 2019 — 210 из 670 ранов
проекта попали в каталог 722 геномов статьи, и почти все они —
конкретные породы вроде Yorkshire Terrier, Bernese Mountain Dog).
Village dogs и дикие псовые в каталоге статьи числятся под другими
BioProject'ами (PRJNA232497, PRJNA233638, PRJNA192935, PRJEB2162).

Из полного списка (130 подходящих геномов после фильтрации) отобран
пилот: 6 village dogs (Vietnam, China ×4) + 6 wild canids (4 Grey Wolf,
Coyote, Dingo), по одному WGS-рану на образец, отобранных по
наименьшему реальному объёму данных (Gbp), а не по сжатому размеру
`.sra` — эти два показателя слабо коррелируют в этом датасете.

**Итог: 11/12 успешно обработаны.** Dingo (SRR2149863) исключён после
6 разных неудачных попыток извлечения (баги fasterq-dump на этом
конкретном файле — `split 3` конфликт со stdout, сегфолт при
`--concatenate-reads`, нехватка места на несжатом fastq; затем
нестабильная сеть при прямой закачке готового fastq.gz с ENA).

Финальная выборка: **6 village dogs + 5 wild canids** (Grey Wolf ×4,
Coyote ×1). Список: `metadata/pilot11_final.tsv`.

## 2. Референс и регионы

- Референс: **CanFam3.1** (GCF_000002285.3), полный геном, bwa-индекс.
- Регионы: **246 CDR** из Table S5 (не Table S4, как значилось в
  исходной формулировке задачи — по факту оказалось, что 246 XP-CLR
  окон лежат в Table S5, Additional File 1; Table S4/Additional File 5
  — другие данные, outlier F_ST loci и aCGH-графики соответственно).
  Проверено: 246 строк, 10.81 Mb суммарно — точно совпадает с
  заявленным в статье "10.8 Mb". RAI1 и RNPC3 присутствуют, как и
  ожидалось.
- BED-файл (координаты переведены из "chr5" в RefSeq accession'ы через
  GFF3 CanFam3.1): `metadata/pendleton2018_246CDR.bed`

## 3. Пайплайн и дисковые проблемы (важный технический урок)

Тот же принцип, что и в fox/fennec-проектах: `prefetch` →
`fasterq-dump --stdout` → `bwa mem` → `samtools view -F 4` →
`samtools sort`, без сырых fastq на диске.

**Ключевая находка, отличавшая этот пилот от предыдущих:** полногеномное
выравнивание (не region-restricted, как в fox-проекте) даёт BAM
размером в несколько GB на образец. После 4 готовых BAM (~19GB)
+ референс (6.7GB) свободного места стало недостаточно для temp-фазы
fasterq-dump у ЛЮБОГО следующего образца — независимо от его размера
(7 провалов подряд с одной и той же `storage exhausted`). Решение:
**обрезка BAM до 246 CDR сразу после выравнивания**, с удалением
полногеномной версии — сократило размер BAM с ~4-5GB до ~20-25MB
(проверено, не оценка) и полностью убрало эффект накопительного "пола".

Также по пути:
- исправлен баг с ненужной передокачкой `.sra` при retry (скрипт
  безусловно чистил рабочую директорию перед каждой попыткой)
- добавлен двухуровневый автовосстанавливающийся супервизор
  (`02_supervisor.sh`) + внешний вотчдог на cron/launchd
  (`03_watchdog.sh`) — среда несколько раз убивала весь процесс-дерево
  целиком под дисковым давлением, включая сам супервизор; вотчдог
  переживает это, т.к. запускается независимо (снят после завершения
  пилота)

## 4. Вызов вариантов (11 образцов, 246 CDR)

`bcftools mpileup -a AD,DP` + `bcftools call -mv`, совместный вызов,
мягкий фильтр `MQ<40 → LowMQ` (тот же критерий, что в fox/fennec).

- Всего сайтов: **61,271**
- PASS: **60,498 (98.7%)**
- LowMQ: **773 (1.3%)**

Доля LowMQ заметно ниже, чем в fox/fennec-проектах (15-28% там) —
246 CDR представляют собой курируемые уникальные (не повторяющиеся)
регионы генома, в отличие от произвольного паддинга вокруг генов.

## 5. Private/shared/core и village vs wild (n=6+5)

Скрипт: `scripts/02_village_vs_wild_analysis.py`

- Private (1 образец): 30,962 (51.2%)
- Shared (2+): 29,536
- Core (все 11): 376

**519 group-fixed сайтов** (ALT фиксирован во всех village ИЛИ во всех
wild, отсутствует в другой группе) — заметно больше по доле, чем в
fox-проекте, что ожидаемо: эти 246 регионов изначально отобраны в
статье именно как наиболее дифференцированные между village dogs и
волками по полногеномному скану, так что повторное обнаружение сильной
дифференциации на них — не открытие, а ожидаемая проверка, что наши
(гораздо более скромные, n=6+5) данные воспроизводят направление
сигнала оригинальной статьи.

### 5.1 Проверка в именованных регионах

**RAI1** (NC_006587.3:41,790,001–41,835,000): 3 group-fixed сайта, все
wild-fixed (41,800,483 / 41,810,535 / 41,824,661).

**RNPC3/AMY2B** (NC_006588.3:47,060,001–47,095,000): 5 group-fixed
сайтов — 3 wild-fixed (47,086,083 / 47,087,134 / 47,088,332), 2
village-fixed (47,086,804 / 47,090,936). AMY2B (амилаза, классическая
находка про адаптацию к крахмалистой диете при доместикации) сама по
себе в основном известна по вариации числа копий (CNV), не SNP —
Pendleton et al. разбирают это отдельным анализом (Additional files
3-5), не в этой XP-CLR таблице; village-fixed SNP здесь могут быть
структурно связаны с этим локусом, но это не прямая проверка CNV.

**Как и везде в этом проекте: n=6+5 — pilot-масштаб, hypothesis-
generating, не публикационная когорта.**

## Файлы

```
village_dogs_data/
├── environment.yaml
├── metadata/
│   ├── runinfo.csv                        — все 670 ранов PRJNA448733
│   ├── selected_village_wild.tsv          — 130 отфильтрованных геномов (весь пул)
│   ├── selected_runinfo.csv               — SRA runinfo для этих 130
│   ├── pilot6_samples.tsv                 — исходный пилот 6+6 (включая Dingo)
│   ├── pilot11_final.tsv                  — финальные 11 (без Dingo)
│   └── pendleton2018_246CDR.bed           — 246 CDR регионов, RefSeq coords
├── supplementary/
│   ├── Supplementary_Data_1.xlsx/.tsv     — Plassais et al. 2019, каталог 722 геномов
│   └── pendleton2018/AdditionalFile1.xlsx — Pendleton et al. 2018, Tables S1-S13
├── reference/                              — CanFam3.1 + bwa-индекс (gitignored-эквивалент, локально)
├── scripts/
│   ├── 01_download_and_align.sh           — скачивание+выравнивание+обрезка до 246 CDR
│   ├── 02_supervisor.sh                   — авто-рестарт при убийстве процесса
│   ├── 02_village_vs_wild_analysis.py     — private/shared/core анализ
│   └── 03_watchdog.sh                     — внешний cron/launchd вотчдог
├── aligned/*.sorted.bam(.bai)              — 11 обрезанных до 246 CDR BAM (~20-25MB каждый)
└── variants/
    ├── village_wild_candidate_regions.vcf              — все сайты, FILTER=PASS/LowMQ
    ├── village_wild_candidate_regions.unfiltered.vcf
    └── village_wild_candidate_regions.PASS.vcf
```
