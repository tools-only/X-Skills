# Index

| Property | Value |
|----------|-------|
| **Name** | Index |
| **Repository** | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/release_notes/v1.65.4-stable/index.md) (🔥 35.1k) |
| **Original Path** | `docs/my-website/release_notes/v1.65.4-stable/index.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | ]
hide_table_of_contents: false
---

import Image from '@theme/IdealImage';
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## Deploy this version

<Tabs>
<TabItem value="docker" label="Docker">

``` showLineNumbers title="docker run litellm"
docker run
-e STORE_MODEL_IN_DB=True
-p 4000:4000
docker.litellm.ai/berriai/litellm:main-v1.65.4-stable
```
</TabItem>

<TabItem value="pip" label="Pip">

``` showLineNumbers title="pip install litellm"
pip install litellm==1.65.4.post1
```
</TabItem>
</Tabs>

v1.65.4-stable is live. Here are the improvements since v1.65.0-stable.

## Key Highlights
- **Preventing DB Deadlocks**: Fixes a high-traffic issue when multiple instances were writing to the DB at the same time. 
- **New Usage Tab**: Enables viewing spend by model and customizing date range

Let's dive in. 

### Preventing DB Deadlocks

<Image img={require('../../img/prevent_deadlocks.jpg')} />

This release fixes the DB deadlocking issue that users faced in high traffic (10K+ RPS). This is great because it enables user/key/team spend tracking works at that scale.

Read more about the new architecture [here |
| **Created** | 2025-04-05 |
| **Updated** | 2025-12-16 |
| **File Hash** | `11b8d378dff64dfe...` |

## Description

tags: []
hide_table_of_contents: false

**Tags:** `]
hide_table_of_contents: false
---

import Image from '@theme/IdealImage';
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## Deploy this version

<Tabs>
<TabItem value="docker" label="Docker">

``` showLineNumbers title="docker run litellm"
docker run
-e STORE_MODEL_IN_DB=True
-p 4000:4000
docker.litellm.ai/berriai/litellm:main-v1.65.4-stable
```
</TabItem>

<TabItem value="pip" label="Pip">

``` showLineNumbers title="pip install litellm"
pip install litellm==1.65.4.post1
```
</TabItem>
</Tabs>

v1.65.4-stable is live. Here are the improvements since v1.65.0-stable.

## Key Highlights
- **Preventing DB Deadlocks**: Fixes a high-traffic issue when multiple instances were writing to the DB at the same time. 
- **New Usage Tab**: Enables viewing spend by model and customizing date range

Let's dive in. 

### Preventing DB Deadlocks

<Image img={require('../../img/prevent_deadlocks.jpg')} />

This release fixes the DB deadlocking issue that users faced in high traffic (10K+ RPS). This is great because it enables user/key/team spend tracking works at that scale.

Read more about the new architecture [here`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/release_notes/v1.65.4-stable/index.md)*
