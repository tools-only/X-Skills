---
name: guidewire-common-errors
description: |
  Diagnose and resolve common Guidewire InsuranceSuite errors.
  Use when encountering API errors, Gosu exceptions, validation failures,
  or integration issues in PolicyCenter, ClaimCenter, or BillingCenter.
  Trigger with phrases like "guidewire error", "policycenter error",
  "claimcenter error", "gosu exception", "api error 422".
allowed-tools: Read, Write, Edit, Bash(curl:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Guidewire Common Errors

## Overview

Diagnose and resolve the most common errors encountered in Guidewire InsuranceSuite Cloud API, Gosu development, and integrations.

## Prerequisites

- Access to Guidewire Cloud Console logs
- Understanding of HTTP status codes
- Familiarity with Gosu exception handling

## HTTP Status Code Reference

| Code | Meaning | Common Cause | Resolution |
|------|---------|--------------|------------|
| 400 | Bad Request | Malformed JSON, invalid field | Check request body structure |
| 401 | Unauthorized | Invalid/expired token | Refresh OAuth token |
| 403 | Forbidden | Missing API role | Check role assignments in GCC |
| 404 | Not Found | Invalid endpoint or ID | Verify URL and resource ID |
| 409 | Conflict | Concurrent modification | Retry with updated checksum |
| 422 | Unprocessable | Business rule violation | Fix validation errors |
| 429 | Too Many Requests | Rate limit exceeded | Implement backoff |
| 500 | Server Error | Internal error | Check server logs |
| 503 | Service Unavailable | Maintenance/overload | Retry with backoff |

## Cloud API Errors

### 401 Authentication Errors

```typescript
// Error: "invalid_token" or "expired_token"
// Solution: Implement token refresh
class TokenManager {
  private token: string | null = null;
  private expiry: Date | null = null;

  async getValidToken(): Promise<string> {
    if (this.isTokenValid()) {
      return this.token!;
    }

    try {
      const response = await this.refreshToken();
      this.token = response.access_token;
      // Set expiry 60 seconds before actual expiry for safety
      this.expiry = new Date(Date.now() + (response.expires_in - 60) * 1000);
      return this.token;
    } catch (error) {
      console.error('Token refresh failed:', error);
      throw new AuthenticationError('Unable to obtain valid token');
    }
  }

  private isTokenValid(): boolean {
    return this.token !== null &&
           this.expiry !== null &&
           this.expiry > new Date();
  }
}
```

### 403 Permission Errors

```json
// Error response
{
  "status": 403,
  "error": "Forbidden",
  "message": "Caller does not have permission to access endpoint: POST /account/v1/accounts",
  "details": {
    "requiredRole": "pc_account_admin",
    "callerRoles": ["pc_policy_read"]
  }
}
```

**Resolution:**
1. Log into Guidewire Cloud Console
2. Navigate to Identity & Access > API Roles
3. Find your service account
4. Add the required API role (`pc_account_admin`)
5. Wait 5 minutes for role propagation

### 409 Conflict (Optimistic Locking)

```typescript
// Error: Checksum mismatch indicates concurrent modification
async function updateWithRetry<T>(
  path: string,
  getData: () => Promise<T>,
  updateFn: (data: T) => Partial<T>,
  maxRetries: number = 3
): Promise<T> {
  let attempt = 0;

  while (attempt < maxRetries) {
    try {
      // Get current state with checksum
      const current = await client.request<{ data: T; checksum: string }>('GET', path);

      // Apply update
      const updated = {
        data: {
          attributes: updateFn(current.data),
          checksum: current.checksum
        }
      };

      return await client.request<T>('PATCH', path, updated);
    } catch (error) {
      if (error.response?.status === 409 && attempt < maxRetries - 1) {
        attempt++;
        console.log(`Conflict detected, retrying (${attempt}/${maxRetries})`);
        await sleep(1000 * attempt); // Exponential backoff
        continue;
      }
      throw error;
    }
  }

  throw new Error('Max retries exceeded for optimistic locking conflict');
}
```

### 422 Validation Errors

```typescript
// Parse and display validation errors
interface ValidationErrorResponse {
  status: 422;
  error: 'Unprocessable Entity';
  details: Array<{
    field: string;
    message: string;
    code: string;
    rejectedValue?: any;
  }>;
}

function handleValidationError(error: ValidationErrorResponse): never {
  console.error('Validation Errors:');

  error.details.forEach(detail => {
    console.error(`  [${detail.code}] ${detail.field}: ${detail.message}`);
    if (detail.rejectedValue !== undefined) {
      console.error(`    Rejected value: ${JSON.stringify(detail.rejectedValue)}`);
    }
  });

  // Group by field for form display
  const fieldErrors = error.details.reduce((acc, detail) => {
    acc[detail.field] = acc[detail.field] || [];
    acc[detail.field].push(detail.message);
    return acc;
  }, {} as Record<string, string[]>);

  throw new ValidationError('Request validation failed', fieldErrors);
}
```

**Common Validation Error Codes:**

| Code | Description | Example |
|------|-------------|---------|
| `required` | Missing required field | `accountHolder is required` |
| `invalid_format` | Wrong data format | `dateOfBirth must be ISO date` |
| `invalid_value` | Value not in allowed set | `state.code 'XX' not valid` |
| `range_exceeded` | Value outside range | `premium exceeds policy limit` |
| `duplicate` | Duplicate entity | `Account number already exists` |
| `reference_not_found` | Invalid foreign key | `producerCode.id not found` |

## Gosu Exceptions

### NullPointerException

```gosu
// Bad: Direct property access
var account = getAccount()
print(account.AccountNumber)  // NPE if account is null

// Good: Null-safe access
var account = getAccount()
if (account != null) {
  print(account.AccountNumber)
}

// Better: Elvis operator
var accountNumber = account?.AccountNumber ?: "Unknown"
```

### QueryException

```gosu
// Bad: Returns multiple rows unexpectedly
var account = Query.make(Account)
  .compare(Account#AccountStatus, Equals, AccountStatus.TC_ACTIVE)
  .select()
  .AtMostOneRow  // Throws if multiple results!

// Good: Handle multiple results explicitly
var accounts = Query.make(Account)
  .compare(Account#AccountStatus, Equals, AccountStatus.TC_ACTIVE)
  .select()
  .toList()

if (accounts.Count > 1) {
  LOG.warn("Multiple active accounts found, using first")
}
var account = accounts.first()
```

### TransactionException

```gosu
// Error: Entity from different bundle
function badUpdate(policy : Policy) {
  Transaction.runWithNewBundle(\bundle -> {
    // ERROR: policy is from outside bundle
    policy.Status = PolicyStatus.TC_INFORCE  // TransactionException!
  })
}

// Good: Add entity to bundle first
function goodUpdate(policy : Policy) {
  Transaction.runWithNewBundle(\bundle -> {
    var p = bundle.add(policy)  // Add to new bundle
    p.Status = PolicyStatus.TC_INFORCE
  })
}
```

### ValidationException

```gosu
// Handle validation during commit
try {
  Transaction.runWithNewBundle(\bundle -> {
    var submission = createSubmission(bundle)
    bundle.commit()
  })
} catch (e : gw.api.database.ValidationException) {
  LOG.error("Validation failed: ${e.Message}")
  e.ValidationResult.Errors.each(\err -> {
    LOG.error("  ${err.FieldPath}: ${err.Message}")
  })
  throw e
}
```

## Integration Errors

### REST API Client Errors

```gosu
// Handle external service errors
uses gw.api.rest.FaultToleranceException
uses gw.api.rest.CircuitBreakerOpenException

class IntegrationService {
  function callExternalService(request : Request) : Response {
    try {
      return _client.post("/api/resource", request)
    } catch (e : CircuitBreakerOpenException) {
      LOG.warn("Circuit breaker open, using fallback")
      return getFallbackResponse()
    } catch (e : FaultToleranceException) {
      if (e.Cause typeis java.net.SocketTimeoutException) {
        LOG.error("External service timeout")
        throw new IntegrationTimeoutException("Service timed out", e)
      }
      LOG.error("Integration error: ${e.Message}")
      throw new IntegrationException("External service error", e)
    }
  }
}
```

### Message Queue Errors

```gosu
// Handle message delivery failures
uses gw.api.messaging.MessageTransport

class MessageHandler {
  function handleDeliveryFailure(message : Message, error : Exception) {
    var retryCount = message.RetryCount ?: 0

    if (retryCount < 3) {
      // Retry with exponential backoff
      message.RetryCount = retryCount + 1
      message.ScheduledSendTime = Date.Now.addMinutes(Math.pow(2, retryCount) as int)
      message.Status = MessageStatus.TC_PENDING
      LOG.info("Scheduling retry ${retryCount + 1} for message ${message.ID}")
    } else {
      // Move to dead letter queue
      message.Status = MessageStatus.TC_ERROR
      LOG.error("Message ${message.ID} failed after ${retryCount} retries: ${error.Message}")
      notifySupport(message, error)
    }
  }
}
```

## Database Errors

### Deadlock Detection

```gosu
// Handle database deadlocks
uses java.sql.SQLException

function executeWithDeadlockRetry<T>(operation() : T, maxRetries : int = 3) : T {
  var attempt = 0

  while (true) {
    try {
      return operation()
    } catch (e : SQLException) {
      if (isDeadlock(e) && attempt < maxRetries) {
        attempt++
        LOG.warn("Deadlock detected, retry ${attempt}/${maxRetries}")
        Thread.sleep(100 * attempt)  // Brief delay
        continue
      }
      throw e
    }
  }
}

private function isDeadlock(e : SQLException) : boolean {
  // Oracle: ORA-00060
  // SQL Server: 1205
  // PostgreSQL: 40P01
  return e.ErrorCode == 60 ||
         e.ErrorCode == 1205 ||
         e.SQLState == "40P01"
}
```

## Error Response Wrapper

```typescript
// Standardized error handling wrapper
interface GuidewireError {
  status: number;
  code: string;
  message: string;
  details?: any;
  traceId?: string;
}

class GuidewireApiError extends Error {
  constructor(
    public status: number,
    public code: string,
    message: string,
    public details?: any,
    public traceId?: string
  ) {
    super(message);
    this.name = 'GuidewireApiError';
  }

  static fromResponse(response: any): GuidewireApiError {
    return new GuidewireApiError(
      response.status,
      response.error || 'UNKNOWN',
      response.message || 'An error occurred',
      response.details,
      response.headers?.['x-gw-trace-id']
    );
  }

  toUserMessage(): string {
    switch (this.code) {
      case 'invalid_token': return 'Your session has expired. Please log in again.';
      case 'forbidden': return 'You do not have permission to perform this action.';
      case 'not_found': return 'The requested resource was not found.';
      case 'validation_error': return 'Please correct the errors in your submission.';
      default: return 'An unexpected error occurred. Please try again.';
    }
  }
}
```

## Resources

- [Cloud API Error Reference](https://docs.guidewire.com/cloud/pc/202503/cloudapica/)
- [Gosu Exception Handling](https://gosu-lang.github.io/)
- [Integration Troubleshooting](https://docs.guidewire.com/education/)

## Next Steps

For debugging techniques, see `guidewire-debug-bundle`.
