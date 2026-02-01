# Common Pitfalls and Solutions

| Property | Value |
|----------|-------|
| **Name** | Common Pitfalls and Solutions |
| **Repository** | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/skills/spring-boot/spring-boot-saga-pattern/references/09-pitfalls-solutions.md) (⭐ 80) |
| **Original Path** | `skills/spring-boot/spring-boot-saga-pattern/references/09-pitfalls-solutions.md` |
| **Category** | productivity |
| **Subcategory** | time-management |
| **Tags** | productivity |
| **Created** | 2025-10-28 |
| **Updated** | 2025-10-28 |
| **File Hash** | `9559cf9a4b834bcf...` |

## Description

@Bean
public ConsumerFactory<String, Object> consumerFactory() {
    Map<String, Object> config = new HashMap<>();
    config.put(ConsumerConfig.ENABLE_AUTO_COMMIT_CONFIG, false); // Manual commit
    return new DefaultKafkaConsumerFactory<>(config);
}

**Tags:** `productivity`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/skills/spring-boot/spring-boot-saga-pattern/references/09-pitfalls-solutions.md)*
