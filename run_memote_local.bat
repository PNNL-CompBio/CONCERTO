call activate whole-cell-gemm

memote report snapshot models/azotobacter_vinelandii_dj/iDT1278.xml^
  --filename=az_snapshot_report.html^
  --experimental=models/azotobacter_vinelandii_dj/data/experiments.yml

memote report snapshot models/rhodosporidium_toruloides_ifo_08804/Rt_IFO0880.xml^
  --filename=rt_snapshot_report.html^
  --experimental=models/rhodosporidium_toruloides_ifo_08804/data/experiments.yml


memote report snapshot models/synechococcus_elongatus_pcc_7942/iJB785.xml^
 --filename=se_snapshot_report.html^
 --experimental=models/synechococcus_elongatus_pcc_7942/data/experiments.yml

call conda deactivate
