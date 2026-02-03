# SLIME Source Code Quick Reference

This document provides quick access to important source code files, functions, and classes in the SLIME codebase.

## Core Entry Points

### Main Training Scripts
**Path**: `train.py`
**Purpose**: Synchronous training entry point
**Key Function**: `train(args)` - Main training loop

**Path**: `train_async.py`
**Purpose**: Asynchronous training entry point for overlapping rollout and training
**Usage**: Use instead of `train.py` for higher GPU utilization

## Configuration and Arguments

### Main Argument Parser
**Path**: `slime/utils/arguments.py`
**Key Functions**:
- `get_slime_extra_args_provider()`: Returns arg provider for SLIME-specific args
- `reset_arg(parser, name, **kwargs)`: Resets/adds Megatron argument defaults
- Internal arg groups:
  - `add_cluster_arguments()`: Ray cluster configuration
  - `add_data_arguments()`: Dataset and data processing
  - `add_rollout_arguments()`: Rollout/inference configuration
  - `add_training_arguments()`: RL training hyperparameters
  - `add_reward_model_arguments()`: Reward model setup
  - `add_evaluation_arguments()`: Evaluation configuration

**Key Parameters Defined**:
- Cluster: `--actor-num-nodes`, `--actor-num-gpus-per-node`, `--rollout-num-gpus`, `--colocate`
- Data: `--prompt-data`, `--input-key`, `--label-key`, `--metadata-key`, `--apply-chat-template`
- Rollout: `--rollout-batch-size`, `--n-samples-per-prompt`, `--rollout-max-response-len`, `--rollout-temperature`
- Training: `--num-rollout`, `--num-steps-per-rollout`, `--global-batch-size`, `--advantage-estimator`
- Reward: `--rm-type`, `--custom-rm-path`
- Evaluation: `--eval-interval`, `--eval-prompt-data`, `--n-samples-per-eval-prompt`

### SGLang Arguments Integration
**Path**: `slime/backends/sglang_utils/arguments.py`
**Key Functions**:
- `add_sglang_arguments(parser)`: Adds SGLang parameters with `--sglang-` prefix
- `validate_args(args)`: Validates SGLang-specific arguments

**Usage**: SGLang params passed as `--sglang-<param_name>` (e.g., `--sglang-mem-fraction-static`)

## Data Types and Structures

### Sample Type
**Path**: `slime/utils/types.py`
**Class**: `Sample`
**Key Fields**:
- `prompt: str` - Input prompt text
- `response: str` - Generated response text
- `tokens: list[int]` - Concatenated prompt + response token IDs
- `response_length: int` - Length of response in tokens
- `reward: float` - Reward score from RM
- `loss_mask: list[int]` - Binary mask for loss calculation (1=compute, 0=mask)
- `truncated: bool` - Whether response was truncated
- `aborted: bool` - Whether generation was aborted
- `label: str` - Ground truth label (for evaluation)
- `metadata: dict` - Custom metadata (session_id, tool_code, etc.)

**Usage**: Core data structure passed through rollout → reward → training pipeline

### Rollout Function Types
**Path**: `slime/rollout/base_types.py`
**Classes**:
- `RolloutFnTrainOutput`: Return type for training rollout
- `RolloutFnEvalOutput`: Return type for evaluation rollout

## Rollout and Generation

### Main Rollout Function
**Path**: `slime/rollout/sglang_rollout.py`
**Key Function**: `generate_rollout(args, rollout_id, data_buffer, evaluation=False) -> list[list[Sample]]`
**Purpose**: Main rollout orchestration using SGLang
**Key Features**:
- Asynchronous generation with asyncio
- Dynamic sampling support
- Partial rollout caching
- DP rank balancing

**Key Classes**:
- `GenerateState`: Singleton for global generation state
  - `tokenizer`: HF tokenizer instance
  - `processor`: HF processor for VLMs
  - `semaphore`: Concurrency control
  - `sampling_params`: Default sampling configuration
  - `dp_counts`: DP rank load balancing
  - Methods:
    - `reset()`: Reset state for new rollout
    - `submit_generate_tasks(samples)`: Submit async generation tasks
    - `dp_rank_context()`: Context manager for DP rank selection

**Key Functions**:
- `async def generate_and_rm_group(args, samples, sampling_params, evaluation)`: Generate and compute rewards for sample group
- `async def generate_single(args, sample, sampling_params, evaluation)`: Generate single sample
- `async def call_sglang(prompt, sampling_params, ...)`: Call SGLang API

**Custom Generation Function Signature**:
```python
async def generate(args, sample: Sample, sampling_params) -> Sample:
    # Custom generation logic
    # Must set: sample.tokens, sample.response_length, sample.response, sample.truncated
    return sample
```

### Reward Models
**Path**: `slime/rollout/rm_hub/`
**Key Functions**:
- `async_rm(args, sample: Sample, **kwargs) -> float`: Async reward computation for single sample
- `batched_async_rm(args, samples: list[Sample], **kwargs) -> list[float]`: Batched reward computation

**Custom RM Function Signature**:
```python
async def reward_func(args, sample: Sample, **kwargs) -> float:
    # Custom reward logic
    return reward_score
```

### Filters and Metrics
**Path**: `slime/rollout/filter_hub/`
**Key Files**:
- `dynamic_sampling_filters.py`: Built-in filters
  - `check_reward_nonzero_std(args, samples, **kwargs)`: Filter for reward diversity
- `base_types.py`: Filter interface definitions
  - `MetricGatherer`: Base class for metric collection
  - `call_dynamic_filter(filter_fn, args, samples)`: Execute filter function

## Training Backends

### Megatron Backend
**Path**: `slime/backends/megatron/`
**Key Files**:
- `arguments.py`: Megatron-specific arguments
- `actor.py`: Actor model implementation
- `trainer.py`: Training loop implementation

**Integration Points**:
- Uses `PYTHONPATH` to find Megatron-LM installation
- Imports `megatron.training.arguments.parse_args`
- Checkpoint loading via `load_checkpoint()`, saving via `save_checkpoint()`
- Custom hooks:
  - `--custom-megatron-init-path`: Initialization hook
  - `--custom-megatron-before-log-prob-hook-path`: Pre-log-prob hook
  - `--custom-megatron-before-train-step-hook-path`: Pre-train-step hook

### FSDP Backend
**Path**: `slime/backends/fsdp/`
**Key Features**:
- Native HuggingFace format support (no weight conversion)
- PyTorch FSDP2 implementation
- Auto-reads architecture from `AutoConfig`

**Key Parameters**:
- `--train-backend fsdp`: Enable FSDP backend
- `--hf-checkpoint`: HF model path (required, used for both rollout and training)
- `--gradient-checkpointing`: Enable gradient checkpointing
- `--fsdp-cpu-offload`: Offload to CPU
- `--context-parallel-size`: Context parallelism (supported)

## Utilities

### Data Management
**Path**: `slime/utils/data.py`
**Class**: `Dataset`
**Purpose**: Manage JSONL datasets with streaming support

### Evaluation Configuration
**Path**: `slime/utils/eval_config.py`
**Classes**:
- `EvalDatasetConfig`: Configuration for eval datasets
- Functions:
  - `build_eval_dataset_configs(args)`: Build eval configs from args
  - `ensure_dataset_list(data_arg)`: Parse dataset specifications

### Processing Utilities
**Path**: `slime/utils/processing_utils.py`
**Key Functions**:
- `load_tokenizer(checkpoint, trust_remote_code)`: Load HF tokenizer
- `load_processor(checkpoint, trust_remote_code)`: Load HF processor (for VLMs)
- `encode_image_for_rollout_engine(image)`: Encode images for rollout

### HTTP Utilities
**Path**: `slime/utils/http_utils.py`
**Key Functions**:
- `async def get(url, **kwargs)`: Async HTTP GET
- `async def post(url, data, **kwargs)`: Async HTTP POST
**Usage**: Communicate with SGLang servers and router

### Async Utilities
**Path**: `slime/utils/async_utils.py`
**Key Functions**:
- `run(coro)`: Run async coroutine in sync context
**Usage**: Bridge between sync and async code

### Logging
**Path**: `slime/utils/logging_utils.py`
**Key Functions**:
- `configure_logger(name, level)`: Configure logger instance

### Misc Utilities
**Path**: `slime/utils/misc.py`
**Key Items**:
- `SingletonMeta`: Metaclass for singleton pattern
- `load_function(path)`: Dynamically load function from module path

## Conversion Tools

### HF to Megatron Conversion
**Path**: `tools/convert_hf_to_torch_dist.py`
**Usage**:
```bash
PYTHONPATH=/root/Megatron-LM python tools/convert_hf_to_torch_dist.py \
    ${MODEL_ARGS[@]} \
    --hf-checkpoint /path/to/hf_model \
    --save /path/to/torch_dist_output
```

### Megatron to HF Conversion
**Path**: `tools/convert_torch_dist_to_hf.py`
**Usage**:
```bash
PYTHONPATH=/root/Megatron-LM python tools/convert_torch_dist_to_hf.py \
    --input-dir /path/to/torch_dist/iter_xxx/ \
    --output-dir /path/to/hf_output \
    --origin-hf-dir /path/to/original_hf_model
```

### FSDP to HF Conversion
**Path**: `tools/convert_fsdp_to_hf.py`
**Usage**:
```bash
python tools/convert_fsdp_to_hf.py \
    --input-dir /path/to/fsdp_ckpt/iter_xxx \
    --output-dir /path/to/hf_output \
    --origin-hf-dir /path/to/original_hf_model
```

## Model Configuration Files

**Path**: `scripts/models/`
**Files**:
- `qwen3-4B.sh`: Defines `MODEL_ARGS` array for Qwen3-4B
- `glm4-9B.sh`: Defines `MODEL_ARGS` array for GLM4-9B
- `qwen3-30B-A3B.sh`: MOE configuration for Qwen3-30B
- `glm45-355B-A32B.sh`: Large MOE configuration

**Key Parameters in Model Configs**:
- Architecture: `--num-layers`, `--hidden-size`, `--ffn-hidden-size`
- Attention: `--num-attention-heads`, `--num-query-groups`, `--kv-channels`
- Normalization: `--normalization RMSNorm`, `--norm-epsilon`
- Position: `--use-rotary-position-embeddings`, `--rotary-base`
- Vocab: `--vocab-size`, `--disable-bias-linear`
- Activation: `--swiglu`

## Quick Code Navigation

### Finding Reward Model Implementations
1. Check `slime/rollout/rm_hub/` for built-in RMs
2. Use `--rm-type <type>` for built-in, `--custom-rm-path` for custom

### Adding Custom Generation Logic
1. Implement `async def generate(args, sample: Sample, sampling_params) -> Sample`
2. Set via `--custom-generate-function-path module.path:function_name`
3. Reference implementation: `slime/rollout/sglang_rollout.py` (default)
4. Multi-turn example: `examples/search-r1/search_r1_logic.py`

### Adding Custom Filters
1. Implement filter function: `def filter_fn(args, samples: list[Sample], **kwargs) -> bool`
2. Set via `--dynamic-sampling-filter-path module.path:function_name`
3. Built-in filters: `slime/rollout/filter_hub/dynamic_sampling_filters.py`

### Customizing Rollout Function
1. Implement `def generate_rollout(args, rollout_id, data_buffer, evaluation) -> list[list[Sample]]`
2. Set via `--rollout-function-path module.path:function_name`
3. Default implementation: `slime/rollout/sglang_rollout.py:generate_rollout`

### Adding Megatron Hooks
- Init hook: `--custom-megatron-init-path module.path:function_name`
- Pre-log-prob hook: `--custom-megatron-before-log-prob-hook-path module.path:function_name`
- Pre-train-step hook: `--custom-megatron-before-train-step-hook-path module.path:function_name`

## Common Import Patterns

```python
# Load SLIME utilities
from slime.utils.types import Sample
from slime.utils.http_utils import post
from slime.utils.processing_utils import load_tokenizer
from slime.rollout.rm_hub import async_rm

# Custom generation function
async def my_generate(args, sample: Sample, sampling_params) -> Sample:
    tokenizer = load_tokenizer(args.hf_checkpoint, trust_remote_code=True)
    # ... generation logic
    return sample

# Custom reward function
async def my_reward(args, sample: Sample, **kwargs) -> float:
    # ... reward computation
    return reward
```
