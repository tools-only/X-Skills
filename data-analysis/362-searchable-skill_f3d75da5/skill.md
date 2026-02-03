# Terminal Bench 2 Skills Index

Searchable index of 89 skills learned from Terminal Bench 2.0. Skills are available in both `trajectory-only/` and `trajectory-feedback/` variants.

## Security & Cryptography

- **break-filter-js-from-html**: Bypassing HTML/JavaScript sanitization filters in security testing contexts, XSS filter bypasses, parser differentials
- **crack-7z-hash**: Cracking 7z archive password hashes, password recovery from encrypted archives
- **feal-differential-cryptanalysis**: Differential cryptanalysis attacks on FEAL and Feistel network ciphers, round key recovery
- **feal-linear-cryptanalysis**: Linear cryptanalysis attacks on FEAL, key recovery from plaintext-ciphertext pairs
- **filter-js-from-html**: Filtering JavaScript and XSS attack vectors from HTML while preserving formatting
- **fix-code-vulnerability**: Identifying and fixing security vulnerabilities, CVEs, CWEs, injection attacks
- **git-leak-recovery**: Recovering secrets from Git history, cleaning sensitive data from repositories
- **model-extraction-relu-logits**: Extracting weight matrices from black-box ReLU neural networks
- **password-recovery**: Digital forensics, recovering passwords from disk images and deleted files
- **sanitize-git-repo**: Removing API keys, tokens, and credentials from codebases
- **vulnerable-secret**: Extracting secrets from protected/obfuscated binaries through reverse engineering

## Machine Learning & AI

- **bn-fit-modify**: Bayesian Networks with pgmpy, structure learning, parameter fitting, do-calculus
- **caffe-cifar-10**: Building Caffe from source, training CIFAR-10 CNN models
- **count-dataset-tokens**: Counting tokens in HuggingFace datasets, tokenizer usage
- **gpt2-codegolf**: Minimal GPT-2 inference, BPE tokenization, checkpoint parsing
- **hf-model-inference**: HuggingFace model inference with Flask APIs
- **llm-inference-batching-scheduler**: Optimizing LLM inference request batching and scheduling
- **mteb-leaderboard**: Querying ML model leaderboards and embedding benchmarks
- **mteb-retrieve**: Text embedding retrieval with sentence transformers, semantic search
- **pytorch-model-cli**: Extracting PyTorch weights and reimplementing inference in C/C++
- **pytorch-model-recovery**: Recovering model architectures from state dictionaries, TorchScript conversion
- **torch-pipeline-parallelism**: PyTorch pipeline parallelism, AFAB scheduling
- **torch-tensor-parallelism**: Tensor parallelism, ColumnParallelLinear, RowParallelLinear
- **train-fasttext**: Training FastText models with size/accuracy constraints

## Bioinformatics & Science

- **dna-assembly**: Golden Gate assembly, Gibson assembly, primer design for Type IIS enzymes
- **dna-insert**: Site-directed mutagenesis primers, Q5 kit, PCR-based insertion
- **mcmc-sampling-stan**: Bayesian MCMC sampling with RStan, hierarchical models
- **protein-assembly**: Designing fusion protein gBlock sequences from PDB databases
- **raman-fitting**: Raman spectrum peak fitting, Lorentzian/Gaussian fitting, graphene characterization
- **rstan-to-pystan**: Converting R-Stan code to Python-Stan (PyStan 3.x)
- **sam-cell-seg**: SAM-based cell segmentation, MobileSAM, mask-to-polygon conversion

## Systems & Build

- **build-cython-ext**: Building Cython extensions, numpy compatibility issues
- **build-pmars**: Building pMARS (Core War simulator) from source
- **build-pov-ray**: Building POV-Ray raytracer from source, legacy C code
- **compile-compcert**: Building CompCert verified C compiler, OCaml/Coq setup
- **custom-memory-heap-crash**: Debugging custom memory heap crashes, use-after-free, static destruction
- **fix-ocaml-gc**: Debugging OCaml garbage collector, runtime C code, pointer arithmetic
- **install-windows-3.11**: Installing Windows 3.11 in QEMU with VNC/noVNC access
- **qemu-alpine-ssh**: QEMU VMs with Alpine Linux and SSH access
- **qemu-startup**: QEMU VM configuration, serial console, telnet access

## Image & Video Processing

- **chess-best-move**: Analyzing chess positions from images, piece detection, move calculation
- **code-from-image**: Extracting code/pseudocode from images using OCR
- **extract-moves-from-video**: Extracting text commands from video recordings using OCR
- **path-tracing**: Implementing path tracers and ray tracers, image reconstruction
- **path-tracing-reverse**: Reverse engineering graphics rendering programs from binaries
- **video-processing**: Video analysis, frame-level event detection with OpenCV

## Git & Version Control

- **configure-git-webserver**: Git repositories with automatic web deployment via post-receive hooks
- **fix-git**: Recovering lost commits, resolving detached HEAD states
- **git-multibranch**: Multi-branch deployment systems with SSH and web servers

## Data Processing

- **financial-document-processor**: Processing invoices/receipts with OCR, data extraction
- **large-scale-text-editing**: Bulk text transformations with Vim, sed, awk
- **log-summary-date-ranges**: Log file analysis, aggregating counts by date and severity
- **multi-source-data-merger**: Merging data from JSON, CSV, Parquet, XML sources
- **reshard-c4-data**: Data resharding, redistributing files across directories

## Web & Infrastructure

- **kv-store-grpc**: gRPC-based key-value store services in Python
- **mailman**: Mailing list servers with Postfix and Mailman3
- **nginx-request-logging**: Nginx configuration, custom logging, rate limiting
- **openssl-selfsigned-cert**: Self-signed SSL/TLS certificates with OpenSSL
- **pypi-server**: Local PyPI servers for hosting Python packages

## Databases

- **db-wal-recovery**: Recovering data from SQLite WAL files, corrupted databases
- **query-optimize**: SQL query optimization, query plan analysis
- **sqlite-db-truncate**: Recovering data from corrupted/truncated SQLite files
- **sqlite-with-gcov**: Compiling SQLite with gcov code coverage

## Embedded & Low-Level

- **circuit-fibsqrt**: Combinational/sequential logic circuits, gate-level descriptions
- **extract-elf**: Parsing ELF binaries, extracting memory contents
- **make-doom-for-mips**: Cross-compiling C programs for embedded MIPS environments
- **make-mips-interpreter**: Building MIPS interpreters/emulators, ELF loaders

## Mathematics & Optimization

- **adaptive-rejection-sampler**: Adaptive Rejection Sampling algorithms
- **constraints-scheduling**: Constraint-based scheduling, calendar conflict resolution
- **distribution-search**: Finding probability distributions satisfying statistical constraints
- **largest-eigenval**: Optimizing eigenvalue computations, beating numpy/scipy
- **portfolio-optimization**: High-performance portfolio optimization with C extensions

## Python & Async

- **cancel-async-tasks**: Asyncio task cancellation, signal handling, graceful shutdown
- **headless-terminal**: Headless terminal interfaces, pseudo-terminal wrappers
- **modernize-scientific-stack**: Migrating Python 2 scientific code to Python 3

## Languages & Compilers

- **cobol-modernization**: Converting COBOL to modern languages (Python, Java)
- **polyglot-c-py**: Creating polyglot files valid in both Python and C
- **polyglot-rust-c**: Creating polyglot source files for Rust and C/C++
- **prove-plus-comm**: Coq proofs, induction, arithmetic lemmas
- **schemelike-metacircular-eval**: Metacircular evaluators, Scheme interpreters

## Games & Simulations

- **regex-chess**: Generating chess moves using only regex pattern matching
- **tune-mjcf**: Optimizing MuJoCo MJCF files for simulation performance
- **winning-avg-corewars**: Developing CoreWars warriors, Redcode assembly

## Text & Documents

- **gcode-to-text**: Extracting text from GCODE files by analyzing toolpath geometry
- **overfull-hbox**: Fixing LaTeX overfull hbox warnings with synonym replacement
- **sparql-university**: SPARQL queries against RDF/Turtle datasets
- **write-compressor**: Implementing encoders compatible with existing decoders

## Other

- **merge-diff-arc-agi-task**: ARC-AGI pattern recognition with git operations
- **regex-log**: Complex regex for log parsing with validation constraints
