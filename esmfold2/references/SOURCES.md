# Primary sources and qualification limits

Retrieved2026-11-13:
- https://github.com/Biohub/esm — published Python API, reference backend, platform notes.
- https://raw.githubusercontent.com/Biohub/esm/main/pyproject.toml — dependency declarations.
- https://pypi.org/pypi/esm/json — release3.4.3, Python>=3.12.
- https://huggingface.co/biohub/ESMFold2/tree/69869f737beffec5294845ede23db5fc0b4f509e — checkpoint revision; approximately26.32GB of F32 weights, bundled language model, additional CCD data.
- https://huggingface.co/biohub/ESMFold2-Fast — approximately26.14GB weights, not materially smaller.
- https://docs.pytorch.org/docs/2.15/generated/torch.linalg.svd.html — non-default driver selection is CUDA-only.

Do not reproduce outdated patches that supply driver="gesvd" on CPU. Current upstream layers.py selects gesvd on CUDA and None otherwise; source review is not proof every MPS operation succeeds.

The script follows the current upstream API but has not run against installed esm3.4.3 or actual weights. Dependency resolution and runtime compatibility must be tested before release promotion. FP32 is a conservative trial configuration, not a universal validated-device guarantee. No published performance or binding claim is made here.
