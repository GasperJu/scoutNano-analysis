| Field | `data2024_triggEff.py` | `mc_triggEff.py` | `mcQCD_triggEff_loop.py` | Recommended authoritative value | Rationale |
|---|---|---|---|---|---|
| Reference trigger path(s) | `PFScouting_SingleMuon` | `PFScouting_SingleMuon` | `PFScouting_SingleMuon` | `PFScouting_SingleMuon` | Already consistent across all producers. |
| Signal trigger path(s) | `PFScouting_JetHT` | `PFScouting_JetHT` | `PFScouting_JetHT` | `PFScouting_JetHT` | Already consistent across all producers. |
| JetID | Eta-region thresholds with NHF/NEMF/muFrac/multiplicity cuts | Same as data | Same as data | Keep data definition | Already harmonized and used in mainstream data path. |
| MuonID quality thresholds | `nRecoMuonMatchedStations>1`, `nValidPixelHits>0`, `nTrackerLayersWithMeasurement>5` | Same as data | Stricter: `>3`, `>1`, `>7` | Keep data definition | Data+main MC match; avoid tightening acceptance unexpectedly. |
| Muon-cleaning | Subtract muon 4-vector from jets for `ΔR<0.4`; drop jets with `pT<1` after subtraction | Same as data | Same as data | Keep data definition | Full overlap subtraction is already common and physics-motivated. |
| L1 seed baseline accept | `DST.PFScouting_JetHT` and (`L1_HTT280er` OR `L1_SingleJet180`) | Same as data | Accepts additional unprescaled bits (`DoubleJet`, `ETT2000`) | Keep data baseline for default selection | Preserves current data behavior and backward compatibility. |
| L1 prescaled/disabled exclusion | Not applied in baseline event acceptance | Not applied in baseline event acceptance | Explicitly rejects prescaled-only (`HTT200/255`) events | Keep as named alternative hypothesis (`exclusive`) in shared config | Enables later systematic studies without changing current baseline behavior. |
| Geometric categorization | Full VBF / Boosted(AK8+ISR) / Resolved(AK4+ISR) / Rest logic | No categorization (inclusive only) | No categorization (inclusive only) | Keep data categorization logic in data producer; MC remains inclusive in this step | Matches existing outputs while still harmonizing T&P object-level definitions. |

## Notes
- Shared module extraction in this step centralizes the reconciled T&P pieces (`REFERENCE_PATHS`, `SIGNAL_PATHS`, `pass_jet_id`, `pass_muon_id`, `clean_jets`) and introduces `L1_SEED_HYPOTHESES` for future multi-hypothesis filling.
- `# TODO(reconcile):` Confirm whether `exclusive` (prescaled-only exclusion) should become default once MC-vs-data closure is validated with multi-hypothesis comparisons.
