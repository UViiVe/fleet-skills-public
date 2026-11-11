# API contract, checked 2026-11-11

The upstream 2.0.0 SDK exposes `alphagenome.models.v2.dna_client`.
`create(api_key, model_version=None, timeout=None)` uses the current server default
when no model version is supplied. This integration does not select AG1 folds.

`model.supported_conditions(Organism.HOMO_SAPIENS)` returns a DataFrame with
one row per supported reference condition and output-specific imputation labels.
`model.predict_interval(interval, organism=..., requested_outputs=[OutputType.RNA_SEQ],
condition=...)` predicts the reference interval for one condition.

The smoke uses 2**14 bases, the smallest supported interval, not the full 1-Mb
variant example in the getting-started guide. It tests transport and returned data,
not the biological adequacy of that context window. Short context can omit distal
regulation. No alternate prediction or effect score is computed.

Sources: https://www.alphagenomedocs.com/api/models_v2.html and
https://www.alphagenomedocs.com/getting_started.html .
