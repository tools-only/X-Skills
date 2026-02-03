# SLIME Documentation Navigation Guide

## Getting Started

### Quick Start (Priority: High)
**Path**: `docs/en/get_started/quick_start.md`
**Purpose**: Essential first read for setting up SLIME within one hour
**Contents**:
- Environment setup (Docker recommended, hardware support)
- Model and dataset download
- Weight conversion (HF ↔ Megatron)
- Training script execution and parameter overview
- Key features: colocated actor/rollout, dynamic sampling, partial rollout, bf16 training with fp8 inference

### Usage Guide (Priority: High)
**Path**: `docs/en/get_started/usage.md`
**Purpose**: Comprehensive parameter reference and configuration guide
**Contents**:
- Cluster resource allocation (actor/rollout GPUs, colocate mode)
- Training backend selection (Megatron vs FSDP)
- Megatron configuration (model params, parallelism, checkpoints)
- SGLang configuration (HF checkpoint, router setup)
- Data format (JSONL with prompt/label/metadata)
- RL hyperparameters (advantage estimators, loss calculation, TIS)
- Custom rollout/reward functions
- SGLang/Megatron integration details
- FSDP backend quick start

### Customization (Priority: Medium)
**Path**: `docs/en/get_started/customization.md`
**Purpose**: Advanced customization of SLIME workflows
**Use when**: Building custom reward models, generation functions, or filters

### Q&A (Priority: Low)
**Path**: `docs/en/get_started/qa.md`
**Purpose**: Common issues and troubleshooting
**Use when**: Debugging setup or runtime issues

## Advanced Features

### Fault Tolerance
**Path**: `docs/en/advanced/fault-tolerance.md`
**Purpose**: Handling failures in distributed training
**Use when**: Setting up multi-node training with automatic recovery

### PD Disaggregation
**Path**: `docs/en/advanced/pd-disaggregation.md`
**Purpose**: Prefill-Decode separation for improved efficiency
**Use when**: Optimizing inference performance with separate prefill/decode servers

### Speculative Decoding
**Path**: `docs/en/advanced/speculative-decoding.md`
**Purpose**: Accelerating inference with draft models
**Use when**: Need faster rollout generation

### Architecture Support Beyond Megatron
**Path**: `docs/en/advanced/arch-support-beyond-megatron.md`
**Purpose**: Using non-Megatron architectures
**Use when**: Working with custom model architectures

## Platform Support

### AMD Tutorial
**Path**: `docs/en/platform_support/amd_tutorial.md`
**Purpose**: Running SLIME on AMD GPUs
**Use when**: Using AMD hardware instead of NVIDIA

## Example Configurations

### Single-Machine Examples
1. **Qwen3-4B** (`docs/en/examples/qwen3-4B.md`): 8xH100, basic setup
2. **GLM4-9B** (`docs/en/examples/glm4-9B.md`): 8xH100, includes dynamic sampling
3. **Qwen3-4B Base (OpenHermes)** (`docs/en/examples/qwen3-4b-base-openhermes.md`): FSDP backend example

### Multi-Node Examples
1. **Qwen3-30B-A3B** (`docs/en/examples/qwen3-30B-A3B.md`): Multi-node MOE training
2. **GLM-4.5-355B-A32B** (`docs/en/examples/glm4.5-355B-A32B.md`): 64xH100, large-scale MOE
3. **DeepSeek-R1** (`docs/en/examples/deepseek-r1.md`): 128xH100, massive scale

## Developer Resources

### Debug Guide
**Path**: `docs/en/developer_guide/debug.md`
**Purpose**: Debugging techniques and tools
**Use when**: Troubleshooting training issues

## Blogs and Release Notes

### Introducing SLIME
**Path**: `docs/en/blogs/introducing_slime.md`
**Purpose**: Overview of SLIME framework and design philosophy

### Release v0.1.0
**Path**: `docs/en/blogs/release_v0.1.0.md`
**Purpose**: Initial release notes and features

## Navigation Tips

1. **For first-time users**: Start with `quick_start.md` → `usage.md`
2. **For basic training**: Read `quick_start.md` and example docs (Qwen3-4B or GLM4-9B)
3. **For advanced features**: Check `advanced/` directory after mastering basics
4. **For custom workflows**: Read `customization.md` and `examples/` directory
5. **For multi-node training**: Study `qwen3-30B-A3B.md` or larger examples
6. **For FSDP backend**: See `usage.md` FSDP section and `qwen3-4b-base-openhermes.md`
