# Choosing a condition and extending responsibly

Run `python scripts/alphagenome2.py conditions --output results/conditions
--confirm-noncommercial` as one shell command after exporting your own key.
Inspect the TSV using a spreadsheet or pandas. Search the human-readable
biosample_name; use the exact condition identifier for follow-on analysis.

Blood CD14-positive monocytes use the composite identifier
`UBERON:0000178_CL:0001054`. A condition's RNA-seq can be observed while another
output is imputed. Do not assume all outputs have the same support.

For REF/ALT predictions follow the official getting-started `predict_variant`
example; for gene-oriented effect statistics use `score_variant` and the documented
v2 scorers. Verify alleles and assembly before requesting either operation. These
methods are intentionally outside the helper's fixed smoke API.

For expression conditioning consult
https://www.alphagenomedocs.com/colabs/v2/conditioning_on_single_cell_data.html .
Use up to 256 human cells, raw integer UMI counts and Ensembl gene IDs. Do not
replace missing consent with publicly available example code. No upload is performed
by this skill's helper or tests.
