# Initial reproducibility receipt

Qualification date: 2026-11-11. Python 3.13.12; SDK alphagenome 2.0.0.
A new isolated environment was used, without changing an existing SDK environment.

Authenticated v2 `supported_conditions(HOMO_SAPIENS)` succeeded: 2,799 rows.
This is a metadata-access result, independently distinct from a prediction.
One human 16,384-bp RNA_SEQ reference-interval prediction also succeeded, returning
finite nonempty values of shape (16384, 12). Interval
`chr2:219235509-219251893:.`, GRCh38; condition `UBERON:0000178_CL:0001054`.
No alternate allele, effect score, patient sequence or single-cell data was sent.

The API used its server default; no immutable model release ID was returned or
asserted. Counts and values may change with service releases. Preserve local JSON
receipts and channel TSVs when repeating the workflow. Qualification is not evidence
for predictive accuracy or universal access with other keys. The helper is not a
hard process deadline and does not promise zero internal SDK transport retries.
