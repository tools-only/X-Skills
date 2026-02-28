# Compliance & Safety (EU AI Act Ready - Aug 2026)

> **Legal Disclaimer**: The Lár framework provides architectural patterns to *assist* with compliance. It does not guarantee compliance on its own. You are responsible for the final validation of your system.

Lár is engineered to meet the stringent requirements of the **EU AI Act (2026)** and **FDA 21 CFR Part 11** for High-Risk AI Systems.

Unlike "Black Box" frameworks that obfuscate decision paths, Lár is a "Glass Box" engine designed for forensic auditability.

---

## EU AI Act Alignment

The EU AI Act (fully enforceable August 2026) imposes strict obligations on "High-Risk" systems (e.g., Medical Devices, Employment, Credit Scoring, Critical Infrastructure).

### 1. Article 12: Record-Keeping (Logging)
**Requirement**: Systems must enable "automatic recording of events ('logs') over their lifetime" to ensure traceability.

**Lár Solution**: `State-Diff Ledger`

Every Lár agent automatically produces a `flight_recorder.json` log. This is not a simple debugging print stream; it is a **forensic ledger** containing:

*   **Timestamp**: UTC-aligned execution time.
*   **Input/Output**: The exact prompt sent and the raw completion received.
*   **Model ID**: The specific version of the model used (e.g., `gpt-4-0613`).
*   **State Diff**: The exact variables changed in memory.

```json
// Example Lár Ledger Entry
{
  "step": 4,
  "node": "TriageNode",
  "timestamp": "2026-08-12T10:00:01Z",
  "state_diff": {
    "risk_score": {"old": 0.1, "new": 0.95},
    "routing": {"old": null, "new": "HUMAN_REVIEW"}
  }
}
```

### 2. Article 13: Transparency & Interpretability
**Requirement**: High-risk AI systems must be designed "in such a way that their operation is sufficiently transparent to enable users to interpret the system's output."

**Lár Solution**: "Glass Box" Architecture

*   **No Hidden Prompts**: Lár does not inject "system prompts" behind your back. You own 100% of the prompt.
*   **Explicit Routing**: The logic flow is defined in standard Python code (Nodes and Edges), not in a hidden neural network or a complex "Agent Executor" loop.
*   **Interpretability**: Any Python developer can read `graph.add_edge("Triage", "Human")` and understand the decision path without needing to understand the LLM's internal weights.

### 3. Article 14: Human Oversight
**Requirement**: Systems must be designed so that they can be "effectively overseen by natural persons," including the ability to "interrupt the system" or "override" decisions.

**Lár Solution**: Native `Interrupt` Pattern

Lár treats "Human Intervention" as a first-class citizen in the graph.

*   **Pause & Resume**: You can execute the graph up to a checkpoint (e.g., `before="ExecuteTool"`), inspect the state, and resume.
*   **State Modification**: A human supervisor can manually edit the memory (e.g., correcting a hallucinated argument) before approving the next step.

---

## FDA 21 CFR Part 11 (Electronic Records)

For healthcare and pharmaceutical applications (e.g., Drug Discovery pipelines), Lár supports the key requirements of **Part 11**:

1.  **Validation**: The deterministic nature of the graph means you can run regression tests. Given Input X and Fixed Seed Y, the graph traverses effectively the same path.
2.  **Audit Trails**: The `GraphExecutor` logs are immutable and time-stamped.
3.  **Authority Checks**: Lár's `SecurityNode` pattern allows you to implement permissions (e.g., "Only User A can approve Tool B") directly in the graph logic.

### Cryptographic Audit Trails (HMAC Signing)

To truly comply with enterprise regulations (like HIPAA, SOC2, FDA GxP, SEC/FINRA), an audit log is not enough—you must prove mathematically that the log has not been tampered with.

Lár v1.5.1 introduced **Cryptographic Signatures** for the audit log. By passing an `hmac_secret` (e.g., from AWS KMS or HashiCorp Vault) to the `GraphExecutor`, the engine will sign the final JSON execution trace using HMAC-SHA256. 

```python
from lar import GraphExecutor

# Instantiating the executor with an HMAC secret turns on Cryptographic Auditing
executor = GraphExecutor(
    log_dir="secure_logs", 
    hmac_secret="your_enterprise_secret_key"
)
```

**How to verify (For Auditors):**
If a single character of the payload (like a node output, reasoning string, or token cost) is altered manually after execution, the signature verification will instantly fail. 

We provide a standalone verification script specifically for Compliance Officers to mathematically prove a log's authenticity:

**Step 1:** Locate the generated JSON audit log (e.g., `secure_logs/run_xyz.json`).
**Step 2:** Obtain the enterprise HMAC Secret Key used during the agent's execution.
**Step 3:** Run the verification script from your terminal:
```bash
python examples/compliance/11_verify_audit_log.py secure_logs/run_xyz.json your_enterprise_secret_key
```

**Outcome:** The script will output either `[+] VERIFICATION SUCCESSFUL` (authentic) or `[-] VERIFICATION FAILED` (tampered).

Lár includes four reference implementations to demonstrate this across different industries:
*   [8_hmac_audit_log.py](../examples/compliance/8_hmac_audit_log.py) (Basic usage)
*   [9_high_risk_trading_hmac.py](../examples/compliance/9_high_risk_trading_hmac.py) (Algorithmic Trading & FINRA)
*   [10_pharma_clinical_trials_hmac.py](../examples/compliance/10_pharma_clinical_trials_hmac.py) (Clinical Data Routing & FDA 21 CFR 11)
*   [11_verify_audit_log.py](../examples/compliance/11_verify_audit_log.py) (Standalone Auditor Script)

---

## Risk Mitigation: Dynamic Graphs (Self-Modifying Code)

Lár v1.1 introduces `DynamicNode`, which allows agents to rewrite their execution topology at runtime. While powerful, "Self-Modifying Code" is traditionally a compliance red flag.

**How Lár mitigates this risk:**

### 1. The "Code-as-Event" Principle
In Lár, a topological change is not a hidden internal state. It is an explicit **Event**.
- The `DynamicNode` outputs a JSON `GraphSpec`.
- This JSON spec is **logged physically** in the audit trail before execution.
- **Auditor Verification**: An auditor can replay the exact moment the agent decided to "add a research step" and verify *why* (based on the context).

### 2. Deterministic Topology Validation
The `TopologyValidator` is a **non-AI, deterministic guardrail**.
- **Allowlists**: It enforces that dynamically spawned nodes can ONLY use tools from a pre-approved list. An agent *cannot* invent a "Delete Database" tool if that function isn't in the Python allowlist.
- **Cycle Prevention**: It mathematically proves that the new subgraph is a DAG (Directed Acyclic Graph) or a bounded loop, preventing "Runaway Agent" scenarios.

### 3. Structural Constraints
The modifications are local. A `DynamicNode` can only swap *itself* or its immediate downstream path. It cannot rewrite history or modify upstream nodes, ensuring "Forward-Only" integrity.

---

## Summary for Auditors

| Feature | Lár Implementation | Compliance Value |
| :--- | :--- | :--- |
| **Determinism** | State Machines vs. Loops | Eliminates "Runaway Agent" risk. |
| **Observability** | JSON Flight Recorder | Meets Art. 12 (Recording). |
| **Control** | Standard HIL Patterns | Meets Art. 14 (Oversight). |
| **Privacy** | Local/Air-Gapped Capable | Meets GDPR / Data Sovereignty. |
