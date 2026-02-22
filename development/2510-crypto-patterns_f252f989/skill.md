---
skill_name: crypto-patterns
skill_category: security
description: Cryptographic patterns for encryption, hashing, key management, and TLS
allowed_tools: [Read, Grep, Glob, Edit, Write, WebSearch]
token_estimate: 1000
version: 1.0
last_updated: 2026-01-13
owner: Security Team
status: active
tags: [security, cryptography, encryption, hashing, keys, anti-pattern, validation]
trigger_files: ["**/*crypto*", "**/*encrypt*", "**/*hash*", "**/*auth*", "**/*token*", "**/*jwt*", "**/*secret*"]
trigger_keywords: [encryption, decryption, hash, bcrypt, argon2, jwt, token, secret, key, certificate, tls, ssl, aes, rsa]
quality_keywords: [anti-pattern, vulnerability, secure, validation, best-practice]
---

# Cryptographic Patterns

Secure cryptographic patterns for encryption, hashing, key management, and transport security.

## Purpose

- Guide correct use of cryptographic primitives
- Prevent common crypto implementation mistakes
- Establish secure defaults for key management and storage

---

## Password Hashing

### Recommended Algorithms

| Algorithm | Use Case | Work Factor |
|-----------|----------|-------------|
| **Argon2id** | New implementations (preferred) | Memory: 64MB, Iterations: 3 |
| **bcrypt** | Widely supported, proven | Cost factor: 12+ |
| **scrypt** | Memory-hard alternative | N=2^14, r=8, p=1 |

### Anti-Patterns

| Pattern | Risk | Fix |
|---------|------|-----|
| Plain text storage | Instant compromise | Use adaptive hashing |
| MD5/SHA1 for passwords | Rainbow tables | Use bcrypt/argon2 |
| Low work factor | Brute force | Increase cost/memory |
| Same salt for all | Rainbow tables | Unique salt per password |

### Correct Implementation

```typescript
// GOOD: Argon2id (preferred)
import argon2 from 'argon2';

async function hashPassword(password: string): Promise<string> {
  return argon2.hash(password, {
    type: argon2.argon2id,
    memoryCost: 65536, // 64 MB
    timeCost: 3,
    parallelism: 4
  });
}

async function verifyPassword(hash: string, password: string): Promise<boolean> {
  return argon2.verify(hash, password);
}

// GOOD: bcrypt (widely supported)
import bcrypt from 'bcrypt';

async function hashPassword(password: string): Promise<string> {
  return bcrypt.hash(password, 12); // Cost factor 12
}

async function verifyPassword(hash: string, password: string): Promise<boolean> {
  return bcrypt.compare(password, hash);
}
```

```typescript
// BAD: Never use these for passwords
const hash = crypto.createHash('md5').update(password).digest('hex');
const hash = crypto.createHash('sha256').update(password).digest('hex');
const hash = crypto.createHash('sha1').update(password).digest('hex');
```

---

## Symmetric Encryption

### Recommended Algorithms

| Algorithm | Mode | Use Case |
|-----------|------|----------|
| **AES-256-GCM** | Authenticated | General encryption (preferred) |
| **ChaCha20-Poly1305** | Authenticated | Mobile/embedded |
| **AES-256-CBC** | Non-authenticated | Legacy (add HMAC) |

### Anti-Patterns

| Pattern | Risk | Fix |
|---------|------|-----|
| ECB mode | Pattern exposure | Use GCM or CBC |
| Hardcoded keys | Key compromise | Use key management |
| Reused IV/nonce | Plaintext recovery | Generate random IV |
| No authentication | Ciphertext tampering | Use GCM or add HMAC |

### Correct Implementation

```typescript
import crypto from 'crypto';

// GOOD: AES-256-GCM (authenticated encryption)
function encrypt(plaintext: string, key: Buffer): string {
  const iv = crypto.randomBytes(12); // 96-bit IV for GCM
  const cipher = crypto.createCipheriv('aes-256-gcm', key, iv);

  let encrypted = cipher.update(plaintext, 'utf8', 'hex');
  encrypted += cipher.final('hex');
  const authTag = cipher.getAuthTag();

  // Return: iv + authTag + ciphertext
  return iv.toString('hex') + authTag.toString('hex') + encrypted;
}

function decrypt(ciphertext: string, key: Buffer): string {
  const iv = Buffer.from(ciphertext.slice(0, 24), 'hex');
  const authTag = Buffer.from(ciphertext.slice(24, 56), 'hex');
  const encrypted = ciphertext.slice(56);

  const decipher = crypto.createDecipheriv('aes-256-gcm', key, iv);
  decipher.setAuthTag(authTag);

  let decrypted = decipher.update(encrypted, 'hex', 'utf8');
  decrypted += decipher.final('utf8');
  return decrypted;
}
```

```typescript
// BAD: ECB mode (patterns visible)
const cipher = crypto.createCipheriv('aes-256-ecb', key, null);

// BAD: Reused IV
const iv = Buffer.alloc(16, 0); // Static IV
const cipher = crypto.createCipheriv('aes-256-cbc', key, iv);

// BAD: No authentication (CBC without HMAC)
const cipher = crypto.createCipheriv('aes-256-cbc', key, iv);
// Missing: HMAC verification of ciphertext
```

---

## Key Management

### Key Generation

```typescript
// GOOD: Cryptographically secure random key
const key = crypto.randomBytes(32); // 256 bits

// GOOD: Key derivation from password
import crypto from 'crypto';

function deriveKey(password: string, salt: Buffer): Buffer {
  return crypto.pbkdf2Sync(password, salt, 100000, 32, 'sha256');
}

// Or use scrypt
function deriveKeyScrypt(password: string, salt: Buffer): Buffer {
  return crypto.scryptSync(password, salt, 32, {
    N: 16384, // CPU/memory cost
    r: 8,     // Block size
    p: 1      // Parallelization
  });
}
```

### Key Storage

| Environment | Storage Method | Access |
|-------------|----------------|--------|
| Development | Environment variables | Process only |
| Production | Secrets manager | IAM-controlled |
| CI/CD | Encrypted secrets | Pipeline only |

```typescript
// GOOD: Load from environment
const key = Buffer.from(process.env.ENCRYPTION_KEY!, 'base64');

// GOOD: AWS Secrets Manager
import { SecretsManager } from '@aws-sdk/client-secrets-manager';
const client = new SecretsManager({});
const secret = await client.getSecretValue({ SecretId: 'app/encryption-key' });
const key = Buffer.from(secret.SecretString!, 'base64');

// BAD: Hardcoded key
const key = Buffer.from('my-secret-key-12345678901234567890');

// BAD: Key in source control
// .env committed to git with ENCRYPTION_KEY=...
```

### Key Rotation

```typescript
// Support multiple key versions
interface EncryptedData {
  keyVersion: number;
  ciphertext: string;
}

async function decrypt(data: EncryptedData): Promise<string> {
  const key = await getKey(data.keyVersion);
  return decryptWithKey(data.ciphertext, key);
}

// Re-encrypt with new key during rotation
async function rotateEncryption(data: EncryptedData, newVersion: number): Promise<EncryptedData> {
  const plaintext = await decrypt(data);
  const newKey = await getKey(newVersion);
  return {
    keyVersion: newVersion,
    ciphertext: encryptWithKey(plaintext, newKey)
  };
}
```

---

## JWT Security

### Signing Algorithms

| Algorithm | Type | Use Case |
|-----------|------|----------|
| **RS256** | Asymmetric | Distributed systems (preferred) |
| **ES256** | Asymmetric | Smaller tokens |
| **HS256** | Symmetric | Single service only |

### Anti-Patterns

| Pattern | Risk | Fix |
|---------|------|-----|
| `alg: none` | No verification | Reject unsigned tokens |
| HS256 with public key | Key confusion | Verify algorithm match |
| No expiration | Token reuse | Set short `exp` |
| Sensitive data in payload | Information leak | Encrypt or omit |

### Correct Implementation

```typescript
import jwt from 'jsonwebtoken';

// GOOD: RS256 with proper validation
const publicKey = fs.readFileSync('public.pem');
const privateKey = fs.readFileSync('private.pem');

function createToken(payload: object): string {
  return jwt.sign(payload, privateKey, {
    algorithm: 'RS256',
    expiresIn: '15m',
    issuer: 'my-app',
    audience: 'my-app-users'
  });
}

function verifyToken(token: string): object {
  return jwt.verify(token, publicKey, {
    algorithms: ['RS256'], // Explicitly allow only RS256
    issuer: 'my-app',
    audience: 'my-app-users'
  });
}
```

```typescript
// BAD: Algorithm confusion vulnerability
jwt.verify(token, key); // No algorithm restriction

// BAD: No expiration
jwt.sign(payload, key); // Token valid forever

// BAD: Sensitive data in payload
jwt.sign({ ssn: '123-45-6789', creditCard: '...' }, key);
```

---

## Random Number Generation

```typescript
// GOOD: Cryptographically secure random
import crypto from 'crypto';

// Random bytes
const bytes = crypto.randomBytes(32);

// Random integer in range [0, max)
function secureRandomInt(max: number): number {
  const range = max;
  const bytesNeeded = Math.ceil(Math.log2(range) / 8);
  let value: number;
  do {
    value = crypto.randomBytes(bytesNeeded).readUIntBE(0, bytesNeeded);
  } while (value >= Math.floor(256 ** bytesNeeded / range) * range);
  return value % range;
}

// Random UUID
const uuid = crypto.randomUUID();
```

```typescript
// BAD: Insecure random (predictable)
const insecure = Math.random(); // NOT cryptographically secure
const token = Math.random().toString(36).substring(7); // Predictable
```

---

## TLS Configuration

### Server Configuration

```typescript
import https from 'https';
import fs from 'fs';

const server = https.createServer({
  key: fs.readFileSync('server.key'),
  cert: fs.readFileSync('server.crt'),

  // Modern TLS configuration
  minVersion: 'TLSv1.2',

  // Prefer server cipher order
  honorCipherOrder: true,

  // Strong cipher suites only
  ciphers: [
    'ECDHE-ECDSA-AES128-GCM-SHA256',
    'ECDHE-RSA-AES128-GCM-SHA256',
    'ECDHE-ECDSA-AES256-GCM-SHA384',
    'ECDHE-RSA-AES256-GCM-SHA384'
  ].join(':')
});
```

### Client Configuration

```typescript
import https from 'https';

// GOOD: Verify certificates
const response = await fetch('https://api.example.com', {
  agent: new https.Agent({
    rejectUnauthorized: true // Default, verify certs
  })
});

// BAD: Disable certificate verification
const agent = new https.Agent({
  rejectUnauthorized: false // NEVER in production
});
```

---

## Hashing (Non-Password)

### Use Cases

| Use Case | Algorithm | Example |
|----------|-----------|---------|
| File integrity | SHA-256 | Checksums |
| Data deduplication | SHA-256 | Content addressing |
| HMAC authentication | HMAC-SHA-256 | API signatures |
| Short identifiers | SHA-256 truncated | Short URLs |

```typescript
import crypto from 'crypto';

// File/data integrity
function hashData(data: Buffer): string {
  return crypto.createHash('sha256').update(data).digest('hex');
}

// HMAC for API authentication
function createHmac(message: string, secret: Buffer): string {
  return crypto.createHmac('sha256', secret).update(message).digest('hex');
}

function verifyHmac(message: string, signature: string, secret: Buffer): boolean {
  const expected = createHmac(message, secret);
  return crypto.timingSafeEqual(
    Buffer.from(signature, 'hex'),
    Buffer.from(expected, 'hex')
  );
}
```

---

## Validation Checklist

### Password Storage
- [ ] Using Argon2id or bcrypt (not SHA/MD5)
- [ ] Cost factor appropriate for hardware
- [ ] Unique salt per password (automatic with bcrypt/argon2)

### Encryption
- [ ] Using authenticated encryption (GCM or ChaCha20-Poly1305)
- [ ] Random IV/nonce per encryption
- [ ] Keys from secure source (not hardcoded)
- [ ] Key rotation plan in place

### JWT/Tokens
- [ ] Algorithm explicitly specified and verified
- [ ] Short expiration times
- [ ] No sensitive data in payload
- [ ] Secure secret/key storage

### TLS
- [ ] TLS 1.2+ required
- [ ] Certificate validation enabled
- [ ] Strong cipher suites only

### Random Generation
- [ ] Using crypto.randomBytes (not Math.random)
- [ ] Sufficient entropy for use case

---

## Related Resources

- [OWASP Cryptographic Storage Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Cryptographic_Storage_Cheat_Sheet.html)
- [OWASP Password Storage Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Password_Storage_Cheat_Sheet.html)
- [Mozilla Server Side TLS](https://wiki.mozilla.org/Security/Server_Side_TLS)
- Related skills: `skill_get("web-security")`
