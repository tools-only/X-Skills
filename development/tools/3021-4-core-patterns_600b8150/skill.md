# The Pattern Library



#### 1. Basic Primitives (`examples/basic/`)
| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`1_simple_triage.py`](../../examples/basic/1_simple_triage.py)** | Classification & Linear Routing |
| **2** | **[`2_reward_code_agent.py`](../../examples/basic/2_reward_code_agent.py)** | Code-First Agent Logic |
| **3** | **[`3_support_helper_agent.py`](../../examples/basic/3_support_helper_agent.py)** | Lightweight Tool Assistant |
| **4** | **[`4_fastapi_server.py`](../../examples/basic/4_fastapi_server.py)** | FastAPI Wrapper (Deploy Anywhere) |

#### 2. Core Patterns (`examples/patterns/`)
| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`1_rag_researcher.py`](../../examples/patterns/1_rag_researcher.py)** | RAG (ToolNode) & State Merging |
| **2** | **[`2_self_correction.py`](../../examples/patterns/2_self_correction.py)** | "Judge" Pattern & Error Loops |
| **3** | **[`3_parallel_execution.py`](../../examples/patterns/3_parallel_execution.py)** | Fan-Out / Fan-In Aggregation |
| **4** | **[`4_structured_output.py`](../../examples/patterns/4_structured_output.py)** | Strict JSON Enforcement |
| **5** | **[`5_multi_agent_handoff.py`](../../examples/patterns/5_multi_agent_handoff.py)** | Multi-Agent Collaboration |
| **6** | **[`6_meta_prompt_optimizer.py`](../../examples/patterns/6_meta_prompt_optimizer.py)** | Self-Modifying Agents (Meta-Reasoning) |
| **7** | **[`7_integration_test.py`](../../examples/patterns/7_integration_test.py)** | Integration Builder (CoinCap) |
| **8** | **[`8_ab_tester.py`](../../examples/patterns/8_ab_tester.py)** | A/B Tester (Parallel Prompts) |
| **9** | **[`9_resumable_graph.py`](../../examples/patterns/9_resumable_graph.py)** | Time Traveller (Crash & Resume) |

#### 3. Compliance & Safety (`examples/compliance/`)
| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`1_human_in_the_loop.py`](../../examples/compliance/1_human_in_the_loop.py)** | User Approval & Interrupts |
| **2** | **[`2_security_firewall.py`](../../examples/compliance/2_security_firewall.py)** | Blocking Jailbreaks with Code |
| **3** | **[`3_juried_layer.py`](../../examples/compliance/3_juried_layer.py)** | Proposer -> Jury -> Kernel |
| **4** | **[`4_access_control_agent.py`](../../examples/compliance/4_access_control_agent.py)** | **Flagship Access Control** |
| **5** | **[`5_context_contamination_test.py`](../../examples/compliance/5_context_contamination_test.py)** | Red Teaming: Social Engineering |
| **6** | **[`6_zombie_action_test.py`](../../examples/compliance/6_zombie_action_test.py)** | Red Teaming: Stale Authority |
| **7** | **[`7_hitl_agent.py`](../../examples/compliance/7_hitl_agent.py)** | Article 14 Compliance Node |
| **8** | **[`8_hmac_audit_log.py`](../../examples/compliance/8_hmac_audit_log.py)** | Cryptographic Audits (HMAC) |
| **9** | **[`9_high_risk_trading_hmac.py`](../../examples/compliance/9_high_risk_trading_hmac.py)** | Financial Ledger Signing |
| **10** | **[`10_pharma_clinical_trials_hmac.py`](../../examples/compliance/10_pharma_clinical_trials_hmac.py)** | FDA/GxP Secured Logs |

#### 4. High Scale (`examples/scale/`)
| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`1_corporate_swarm.py`](../../examples/scale/1_corporate_swarm.py)** | **Stress Test**: 60+ Node Graph |
| **2** | **[`2_mini_swarm_pruner.py`](../../examples/scale/2_mini_swarm_pruner.py)** | Dynamic Graph Pruning |
| **3** | **[`3_parallel_newsroom.py`](../../examples/scale/3_parallel_newsroom.py)** | True Parallelism (`BatchNode`) |
| **4** | **[`4_parallel_corporate_swarm.py`](../../examples/scale/4_parallel_corporate_swarm.py)** | Concurrent Branch Execution |

#### 5. Metacognition (`examples/metacognition/`)

See the **[Metacognition Docs](9-metacognition.md)** for a deep dive.

| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`1_dynamic_depth.py`](../../examples/metacognition/1_dynamic_depth.py)** | **Adaptive Complexity** (1 Node vs N Nodes) |
| **2** | **[`2_tool_inventor.py`](../../examples/metacognition/2_tool_inventor.py)** | **Self-Coding** (Writing Tools at Runtime) |
| **3** | **[`3_self_healing.py`](../../examples/metacognition/3_self_healing.py)** | **Error Recovery** (Injecting Fix Subgraphs) |
| **4** | **[`4_adaptive_deep_dive.py`](../../examples/metacognition/4_adaptive_deep_dive.py)** | **Recursive Research** (Spawning Sub-Agents) |
| **5** | **[`5_expert_summoner.py`](../../examples/metacognition/5_expert_summoner.py)** | **Dynamic Persona Instantiation** |

#### 6. Advanced Showcase (`examples/advanced/`)
| # | Pattern | Concept |
| :---: | :--- | :--- |
| **1** | **[`fractal_polymath.py`](../../examples/advanced/fractal_polymath.py)** | **Fractal Agency** (Recursion + Parallelism) |
