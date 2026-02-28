## Privacy Policy

This privacy policy describes how **Skill Agent** (the “Plugin”) processes data when used in Dify.

### 1. Who We Are

- Author: lfenghx
- Repository: <https://github.com/lfenghx/skill_agent>
- Contact:
  - Email: 550916599@qq.com
  - GitHub: lfenghx
  - Bilibili: 元视界\_O凌枫o

### 2. Data We Process

The Plugin may process the following data to provide its functionality:

- **User input**: the `query` parameter and related conversation context passed by Dify.
- **Model selection/configuration**: the `model` selector and parameters used to call the LLM through Dify’s model runtime.
- **Generated artifacts**: files created in the Plugin’s temporary session directory during execution (e.g., `.txt`, `.md`, `.pdf`, images).
- **Operational logs**: debug logs printed by the Plugin runtime for troubleshooting (may include tool call names, file paths under the temp directory, and execution status).

The Plugin does **not** intentionally collect personal data beyond what is required to execute the user’s request.

### 3. How Data Is Used

Data is used strictly for:

- Selecting and invoking skills from the local `skills/` directory.
- Reading skill documentation (`SKILL.md`) and related skill files as needed.
- Running whitelisted commands (e.g., `python`, `node`) inside controlled directories to generate deliverables.
- Returning the final text and generated files back to Dify as tool outputs.

### 4. Data Sharing & Third Parties

Depending on your Dify configuration, data may be transmitted to:

- **LLM providers configured in Dify**: The Plugin invokes the LLM via Dify’s model runtime. Prompts and context may be sent to the configured provider to generate responses and tool plans.
- **Dify remote install/debug service (optional)**: If you use remote debugging install, plugin installation metadata may be exchanged with the configured remote install server.

The Plugin itself does not add additional third-party analytics/telemetry services.

### 5. Storage & Retention

- **Temporary files**: The Plugin creates a per-run temp session directory under the plugin workspace (e.g., `temp/dify-skill-xxxx/`). These files are used to assemble deliverables. The Plugin may clean up old sessions automatically based on its internal retention logic.
- **Conversation summary state (Dify storage)**: The Plugin may store a compact conversation summary and resume state in Dify’s provided storage to support multi-turn runs.

Retention is primarily controlled by:

- Your deployment environment (filesystem retention/backups)
- Your Dify configuration and storage lifecycle policies

### 6. Security

To reduce security risks:

- The Plugin restricts command execution to a whitelist of executables.
- File reads/writes are constrained to approved directories (skill folder and temp session folder).

However, you should still treat generated files and logs as potentially sensitive if your inputs contain sensitive content.

### 7. Your Choices

- Avoid submitting sensitive personal information or secrets in `query` unless necessary.
- Manage/clear conversation data via Dify if your deployment requires data minimization.
- Remove plugin temp directories from disk if you need immediate cleanup in self-hosted deployments.

### 8. Changes to This Policy

This policy may be updated as the Plugin evolves. Updates will be published in the repository.
