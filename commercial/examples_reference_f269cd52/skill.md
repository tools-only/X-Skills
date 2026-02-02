# SLIME Examples Quick Reference

This document provides quick access to example implementations and scripts in the SLIME repository.

## Core Training Examples (Single-Turn)

### Basic GRPO Training
**Script**: `scripts/run-qwen3-4B.sh`, `scripts/run-glm4-9B.sh`
**Purpose**: Standard GRPO training on math reasoning tasks
**Key Features**:
- Rollout-train loop configuration
- Dynamic batch sizing
- Evaluation setup
- Parameter optimization

**Key Ideas**:
- Use `--rollout-batch-size × --n-samples-per-prompt = --global-batch-size × --num-steps-per-rollout`
- Enable `--use-dynamic-batch-size` with `--max-tokens-per-gpu` for efficiency
- Configure `--advantage-estimator grpo` for GRPO algorithm
- Use `--rm-type deepscaler` for built-in reward model

### Dynamic Sampling (DAPO-style)
**Script**: Example in `docs/en/examples/qwen3-4B.md` and `docs/en/get_started/quick_start.md`
**Purpose**: Filter low-quality samples during rollout
**Key Parameters**:
```bash
--over-sampling-batch-size 64 \
--rollout-batch-size 32 \
--dynamic-sampling-filter-path slime.rollout.filter_hub.dynamic_sampling_filters.check_reward_nonzero_std
```

**Key Ideas**:
- Oversample prompts, filter based on reward diversity
- Automatically trigger new sampling when insufficient valid samples
- Use custom filters via `--dynamic-sampling-filter-path`

### Partial Rollout
**Purpose**: Cache and resume aborted requests for efficiency
**Key Parameters**:
```bash
--partial-rollout \
--buffer-filter-path <custom_filter_path>  # default: pop_first
```

**Key Ideas**:
- Reduces waste during dynamic sampling
- Stores partial generations in buffer for next rollout
- Custom extraction strategies via `--buffer-filter-path`

## Advanced Examples

### 1. DrGRPO Algorithm
**Path**: `examples/DrGRPO/`
**Purpose**: Custom reducer implementation for Dr.GRPO
**Use When**: Implementing custom RL algorithms with specialized reduction logic

### 2. Multi-Turn Interaction (Search-R1)
**Path**: `examples/search-r1/`
**Purpose**: Minimal reproduction of Search-R1 with tool calling
**Key Features**:
- Multi-turn conversation handling
- Tool execution and observation
- Loss masking for tool outputs
- Custom generation function

**Key Ideas**:
- Use `--custom-generate-function-path` for multi-turn logic
- Set `loss_mask=1` for model actions, `loss_mask=0` for tool outputs
- Store tool definitions in `metadata` field
- Implement action parsing and execution loop

**Script Location**: `examples/search-r1/search_r1_logic.py`

### 3. Tau-Bench (Agentic Multi-Turn)
**Path**: `examples/tau-bench/`
**Purpose**: Training in multi-turn tool-use environment
**Key Features**: Similar to Search-R1 but with Tau-bench environment

### 4. ReTool (Tool-Enabled Generation)
**Path**: `examples/retool/`
**Purpose**: Tool-enabled language model generation
**Use When**: Need models to use external tools during generation

### 5. Multi-Agent RL
**Path**: `examples/multi_agent/`
**Purpose**: Running multiple agents in RL training
**Use When**: Implementing multi-agent scenarios

### 6. VLM Training (Vision-Language Models)

#### Single-Turn VLM (GEO3K)
**Path**: `examples/geo3k_vlm/`
**Purpose**: Train VLMs with FSDP on GEO3K dataset
**Key Features**: FSDP backend, single-turn reasoning

#### Multi-Turn VLM (GEO3K)
**Path**: `examples/geo3k_vlm_multi_turn/`
**Purpose**: Multi-turn VLM training on GEO3K
**Key Features**: FSDP backend, multi-turn interactions

#### True On-Policy VLM (Qwen3-VL)
**Path**: `examples/true_on_policy_vlm/`
**Purpose**: Strictly equal log probabilities between inference and training
**Use When**: Need exact on-policy training for VLMs

### 7. Low Precision Training
**Path**: `examples/low_precision/`
**Purpose**: FP8 training and inference
**Key Features**:
- Improved throughput and stability
- bf16 training with fp8 inference
- Download FP8 model variant (e.g., Qwen3-4B-FP8)

**Key Ideas**:
```bash
--hf-checkpoint /path/to/Qwen3-4B-FP8  # Use FP8 variant for rollout
--ref-load /path/to/bf16_torch_dist     # Keep bf16 for training
```

### 8. On-Policy Distillation
**Path**: `examples/on_policy_distillation/`
**Purpose**: Teacher-student distillation within on-policy training
**Use When**: Want to distill a larger teacher model during RL training

### 9. True On-Policy Mode
**Path**: `examples/true_on_policy/`
**Purpose**: Ensure strictly equal log probabilities between inference and training
**Key Parameter**: `--true-on-policy-mode`
**Use When**: Need exact on-policy guarantees

### 10. Fully Async Rollout
**Path**: `examples/fully_async/`
**Purpose**: Fully asynchronous rollout generation
**Key Features**: Higher efficiency through async operations
**Script**: Use `train_async.py` instead of `train.py`

### 11. Formal Math Reasoning
**Path**: `examples/formal_math/`
**Purpose**: Formal math reasoning tasks
**Use When**: Training on formal proof generation or verification

### 12. Multi-Task Evaluation
**Path**: `examples/eval_multi_task/`
**Purpose**: OOD evaluation on multiple tasks (GPQA, IFBench, etc.)
**Use When**: Need to evaluate on diverse benchmarks

### 13. Evaluation with NeMo-Skills
**Path**: `examples/eval/`
**Purpose**: Evaluation environment setup using NeMo-Skills
**Use When**: Need standardized evaluation infrastructure

### 14. Reproducibility
**Path**: `examples/reproducibility/`
**Purpose**: Achieving bitwise experiment reproduction
**Use When**: Need deterministic training for research

### 15. Train-Infer Mismatch Helper
**Path**: `examples/train_infer_mismatch_helper/`
**Purpose**: Algorithmic correction methods (TIS, MIS)
**Use When**: Handling train-inference distribution mismatch

### 16. Strands-SGLang Integration
**Path**: `examples/strands_sglang/`
**Purpose**: Integration with Strands-Agents framework
**Use When**: Using Strands scaffolding for agent development

## Common Training Patterns

### Pattern 1: Colocated Training
```bash
--actor-num-nodes 1 \
--actor-num-gpus-per-node 8 \
--colocate \
--sglang-mem-fraction-static 0.7  # Reduce SGLang memory for Megatron
```

### Pattern 2: Disaggregated Training
```bash
--actor-num-nodes 1 \
--actor-num-gpus-per-node 4 \
--rollout-num-gpus 4 \
--rollout-num-gpus-per-engine 2
```

### Pattern 3: Multi-Node Training
```bash
# Start Ray cluster first
ray start --head --node-ip-address ${MASTER_ADDR} --num-gpus 8
# On other nodes:
ray start --address=${MASTER_ADDR}:6379 --num-gpus 8

# Submit job
ray job submit --address="http://127.0.0.1:8265" \
  --runtime-env-json='{"env_vars": {"PYTHONPATH": "/root/Megatron-LM/"}}' \
  -- python3 train.py --actor-num-nodes 8 --actor-num-gpus-per-node 8 ...
```

### Pattern 4: FSDP Backend
```bash
--train-backend fsdp \
--hf-checkpoint /path/to/hf_model  # No weight conversion needed
```

## Model Configuration Files

Model configs are located in `scripts/models/`:
- `qwen3-4B.sh`: Qwen3-4B architecture params
- `glm4-9B.sh`: GLM4-9B architecture params
- `qwen3-30B-A3B.sh`: Qwen3-30B MOE params
- `glm45-355B-A32B.sh`: GLM-4.5-355B large MOE params

**Usage**: `source scripts/models/qwen3-4B.sh` to load `MODEL_ARGS` array

## Dataset Examples

### Standard Format (JSONL)
```jsonl
{
  "prompt": [{"role": "user", "content": "Question here"}],
  "label": "Expected answer",
  "metadata": {"session_id": "123", "custom_field": "value"}
}
```

### Key Parameters
```bash
--prompt-data /path/to/data.jsonl \
--input-key prompt \
--label-key label \
--metadata-key metadata \
--apply-chat-template  # If prompt uses OpenAI message format
```

## RL Algorithm Selection

Available via `--advantage-estimator`:
- `grpo`: Group Relative Policy Optimization
- `gspo`: Group Soft Policy Optimization
- `reinforce_plus_plus`: Reinforce++ algorithm
- `reinforce_plus_plus_baseline`: Reinforce++ with baseline
- `ppo`: Proximal Policy Optimization
- `on_policy_distillation`: On-policy teacher-student distillation

## Reward Model Types

Built-in types via `--rm-type`:
- `deepscaler`: DeepScaler reward model
- Custom: Use `--custom-rm-path <module.path:function_name>`

## Quick Recipe Finder

1. **Basic math RL training**: `scripts/run-qwen3-4B.sh`
2. **Multi-turn tool use**: `examples/search-r1/` or `examples/tau-bench/`
3. **Vision-language RL**: `examples/geo3k_vlm/`
4. **Large-scale MOE**: `docs/en/examples/glm4.5-355B-A32B.md`
5. **Custom generation logic**: See `examples/search-r1/search_r1_logic.py`
6. **Custom reward function**: Implement `async def reward_func(args, sample: Sample) -> float`
7. **Dynamic sampling**: Add `--over-sampling-batch-size` + `--dynamic-sampling-filter-path`
8. **FSDP instead of Megatron**: Add `--train-backend fsdp`
