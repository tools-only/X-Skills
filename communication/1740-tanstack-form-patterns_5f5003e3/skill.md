# TanStack Form Patterns

Form handling with TanStack Form v1 and Zod validation.

## Basic Form Setup

### With Zod Validation

```tsx
'use client';

import { useForm } from '@tanstack/react-form';
import { zodValidator } from '@tanstack/zod-form-adapter';
import { z } from 'zod';

const userSchema = z.object({
  name: z.string().min(1, 'Name is required').max(100),
  email: z.string().email('Invalid email address'),
  age: z.number().min(18, 'Must be 18 or older').optional(),
});

type UserFormData = z.infer<typeof userSchema>;

export function UserForm({ onSubmit }: { onSubmit: (data: UserFormData) => Promise<void> }) {
  const form = useForm({
    defaultValues: {
      name: '',
      email: '',
      age: undefined,
    } satisfies UserFormData,
    validatorAdapter: zodValidator(),
    validators: {
      onChange: userSchema,
    },
    onSubmit: async ({ value }) => {
      await onSubmit(value);
    },
  });

  return (
    <form
      onSubmit={(e) => {
        e.preventDefault();
        form.handleSubmit();
      }}
    >
      <form.Field name="name">
        {(field) => (
          <div>
            <label htmlFor={field.name}>Name</label>
            <input
              id={field.name}
              value={field.state.value}
              onChange={(e) => field.handleChange(e.target.value)}
              onBlur={field.handleBlur}
            />
            {field.state.meta.isTouched && field.state.meta.errors.length > 0 && (
              <span className="text-red-500 text-sm">{field.state.meta.errors[0]}</span>
            )}
          </div>
        )}
      </form.Field>

      <form.Field name="email">
        {(field) => (
          <div>
            <label htmlFor={field.name}>Email</label>
            <input
              id={field.name}
              type="email"
              value={field.state.value}
              onChange={(e) => field.handleChange(e.target.value)}
              onBlur={field.handleBlur}
            />
            {field.state.meta.isTouched && field.state.meta.errors.length > 0 && (
              <span className="text-red-500 text-sm">{field.state.meta.errors[0]}</span>
            )}
          </div>
        )}
      </form.Field>

      <button
        type="submit"
        disabled={!form.state.canSubmit || form.state.isSubmitting}
      >
        {form.state.isSubmitting ? 'Submitting...' : 'Submit'}
      </button>
    </form>
  );
}
```

---

## Field Array Pattern

### Dynamic List of Items

```tsx
const marketSchema = z.object({
  name: z.string().min(1),
  outcomes: z
    .array(
      z.object({
        label: z.string().min(1, 'Outcome label required'),
        probability: z.number().min(0).max(100),
      })
    )
    .min(2, 'At least 2 outcomes required')
    .refine(
      (outcomes) => {
        const total = outcomes.reduce((sum, o) => sum + o.probability, 0);
        return Math.abs(total - 100) < 0.01;
      },
      { message: 'Probabilities must sum to 100%' }
    ),
});

export function MarketForm() {
  const form = useForm({
    defaultValues: {
      name: '',
      outcomes: [
        { label: 'Yes', probability: 50 },
        { label: 'No', probability: 50 },
      ],
    },
    validatorAdapter: zodValidator(),
    validators: { onChange: marketSchema },
  });

  return (
    <form onSubmit={(e) => { e.preventDefault(); form.handleSubmit(); }}>
      <form.Field name="name">
        {(field) => (
          <input
            value={field.state.value}
            onChange={(e) => field.handleChange(e.target.value)}
            placeholder="Market name"
          />
        )}
      </form.Field>

      <form.Field name="outcomes" mode="array">
        {(field) => (
          <div>
            <h3>Outcomes</h3>
            {field.state.value.map((_, index) => (
              <div key={index} className="flex gap-2">
                <form.Field name={`outcomes[${index}].label`}>
                  {(subField) => (
                    <input
                      value={subField.state.value}
                      onChange={(e) => subField.handleChange(e.target.value)}
                      placeholder="Label"
                    />
                  )}
                </form.Field>

                <form.Field name={`outcomes[${index}].probability`}>
                  {(subField) => (
                    <input
                      type="number"
                      value={subField.state.value}
                      onChange={(e) => subField.handleChange(Number(e.target.value))}
                      placeholder="Probability"
                    />
                  )}
                </form.Field>

                {field.state.value.length > 2 && (
                  <button
                    type="button"
                    onClick={() => field.removeValue(index)}
                  >
                    Remove
                  </button>
                )}
              </div>
            ))}

            <button
              type="button"
              onClick={() => field.pushValue({ label: '', probability: 0 })}
            >
              Add Outcome
            </button>

            {field.state.meta.errors.length > 0 && (
              <span className="text-red-500">{field.state.meta.errors[0]}</span>
            )}
          </div>
        )}
      </form.Field>

      <button type="submit" disabled={!form.state.canSubmit}>
        Create Market
      </button>
    </form>
  );
}
```

---

## Async Validation

### Debounced API Check

```tsx
import { z } from 'zod';

const usernameSchema = z.string()
  .min(3, 'Username must be at least 3 characters')
  .max(20, 'Username must be at most 20 characters')
  .regex(/^[a-zA-Z0-9_]+$/, 'Only letters, numbers, and underscores');

export function UsernameField() {
  const form = useForm({
    defaultValues: { username: '' },
    validatorAdapter: zodValidator(),
  });

  return (
    <form.Field
      name="username"
      validators={{
        onChange: usernameSchema,
        onChangeAsyncDebounceMs: 500,
        onChangeAsync: async ({ value }) => {
          // Skip if basic validation fails
          const syncResult = usernameSchema.safeParse(value);
          if (!syncResult.success) return undefined;

          // Check availability
          const response = await fetch(`/api/check-username?u=${value}`);
          const { available } = await response.json();
          return available ? undefined : 'Username is already taken';
        },
      }}
    >
      {(field) => (
        <div>
          <input
            value={field.state.value}
            onChange={(e) => field.handleChange(e.target.value)}
            onBlur={field.handleBlur}
          />
          {field.state.meta.isValidating && <span>Checking...</span>}
          {field.state.meta.errors.length > 0 && (
            <span className="text-red-500">{field.state.meta.errors[0]}</span>
          )}
        </div>
      )}
    </form.Field>
  );
}
```

---

## Integration with TanStack Query

### Submit via Mutation

```tsx
import { useMutation, useQueryClient } from '@tanstack/react-query';
import { useForm } from '@tanstack/react-form';

export function CreateUserForm() {
  const queryClient = useQueryClient();

  const createMutation = useMutation({
    mutationFn: createUser,
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['users'] });
    },
  });

  const form = useForm({
    defaultValues: { name: '', email: '' },
    onSubmit: async ({ value }) => {
      await createMutation.mutateAsync(value);
    },
  });

  return (
    <form onSubmit={(e) => { e.preventDefault(); form.handleSubmit(); }}>
      {/* Fields */}

      {createMutation.isError && (
        <div className="text-red-500">
          Failed to create user: {createMutation.error.message}
        </div>
      )}

      <button
        type="submit"
        disabled={form.state.isSubmitting || createMutation.isPending}
      >
        {createMutation.isPending ? 'Creating...' : 'Create User'}
      </button>
    </form>
  );
}
```

---

## Multi-Step Form (Wizard)

```tsx
import { useState } from 'react';
import { useForm } from '@tanstack/react-form';

const steps = ['personal', 'contact', 'review'] as const;
type Step = typeof steps[number];

export function MultiStepForm() {
  const [currentStep, setCurrentStep] = useState<Step>('personal');

  const form = useForm({
    defaultValues: {
      firstName: '',
      lastName: '',
      email: '',
      phone: '',
    },
    onSubmit: async ({ value }) => {
      await submitForm(value);
    },
  });

  const goNext = () => {
    const currentIndex = steps.indexOf(currentStep);
    if (currentIndex < steps.length - 1) {
      setCurrentStep(steps[currentIndex + 1]);
    }
  };

  const goBack = () => {
    const currentIndex = steps.indexOf(currentStep);
    if (currentIndex > 0) {
      setCurrentStep(steps[currentIndex - 1]);
    }
  };

  return (
    <form onSubmit={(e) => { e.preventDefault(); form.handleSubmit(); }}>
      {currentStep === 'personal' && (
        <div>
          <h2>Personal Information</h2>
          <form.Field name="firstName">
            {(field) => (
              <input
                value={field.state.value}
                onChange={(e) => field.handleChange(e.target.value)}
                placeholder="First Name"
              />
            )}
          </form.Field>
          <form.Field name="lastName">
            {(field) => (
              <input
                value={field.state.value}
                onChange={(e) => field.handleChange(e.target.value)}
                placeholder="Last Name"
              />
            )}
          </form.Field>
        </div>
      )}

      {currentStep === 'contact' && (
        <div>
          <h2>Contact Information</h2>
          <form.Field name="email">
            {(field) => (
              <input
                type="email"
                value={field.state.value}
                onChange={(e) => field.handleChange(e.target.value)}
                placeholder="Email"
              />
            )}
          </form.Field>
          <form.Field name="phone">
            {(field) => (
              <input
                value={field.state.value}
                onChange={(e) => field.handleChange(e.target.value)}
                placeholder="Phone"
              />
            )}
          </form.Field>
        </div>
      )}

      {currentStep === 'review' && (
        <div>
          <h2>Review</h2>
          <p>Name: {form.state.values.firstName} {form.state.values.lastName}</p>
          <p>Email: {form.state.values.email}</p>
          <p>Phone: {form.state.values.phone}</p>
        </div>
      )}

      <div className="flex gap-2">
        {currentStep !== 'personal' && (
          <button type="button" onClick={goBack}>Back</button>
        )}
        {currentStep !== 'review' ? (
          <button type="button" onClick={goNext}>Next</button>
        ) : (
          <button type="submit" disabled={form.state.isSubmitting}>
            Submit
          </button>
        )}
      </div>
    </form>
  );
}
```

---

## Form State Utilities

### Dirty State Detection

```tsx
function FormWithDirtyCheck() {
  const form = useForm({ defaultValues: { name: '' } });

  const isDirty = form.state.isDirty;

  // Warn on navigation
  useEffect(() => {
    const handleBeforeUnload = (e: BeforeUnloadEvent) => {
      if (isDirty) {
        e.preventDefault();
        e.returnValue = '';
      }
    };
    window.addEventListener('beforeunload', handleBeforeUnload);
    return () => window.removeEventListener('beforeunload', handleBeforeUnload);
  }, [isDirty]);

  return (
    <form>
      {/* fields */}
      <button
        type="button"
        onClick={() => form.reset()}
        disabled={!isDirty}
      >
        Reset
      </button>
    </form>
  );
}
```

### Field-Level Meta

```tsx
<form.Field name="email">
  {(field) => (
    <div>
      <input
        value={field.state.value}
        onChange={(e) => field.handleChange(e.target.value)}
        onBlur={field.handleBlur}
        className={field.state.meta.isTouched && field.state.meta.errors.length > 0 ? 'border-red-500' : ''}
      />

      {/* Show validation status */}
      {field.state.meta.isValidating && <span>Validating...</span>}

      {/* Show error only after touch */}
      {field.state.meta.isTouched && field.state.meta.errors.length > 0 && (
        <span className="text-red-500">{field.state.meta.errors[0]}</span>
      )}

      {/* Show success state */}
      {field.state.meta.isTouched &&
       !field.state.meta.isValidating &&
       field.state.meta.errors.length === 0 && (
        <span className="text-green-500">✓</span>
      )}
    </div>
  )}
</form.Field>
```

---

## Reusable Field Components

```tsx
// components/form-field.tsx
interface FormFieldProps<TFormData> {
  form: FormApi<TFormData>;
  name: keyof TFormData & string;
  label: string;
  type?: 'text' | 'email' | 'password' | 'number';
}

export function FormField<TFormData>({
  form,
  name,
  label,
  type = 'text',
}: FormFieldProps<TFormData>) {
  return (
    <form.Field name={name}>
      {(field) => (
        <div className="space-y-1">
          <label htmlFor={name} className="text-sm font-medium">
            {label}
          </label>
          <input
            id={name}
            type={type}
            value={field.state.value as string}
            onChange={(e) => {
              const value = type === 'number'
                ? Number(e.target.value)
                : e.target.value;
              field.handleChange(value as TFormData[typeof name]);
            }}
            onBlur={field.handleBlur}
            className={`w-full px-3 py-2 border rounded ${
              field.state.meta.isTouched && field.state.meta.errors.length > 0
                ? 'border-red-500'
                : 'border-gray-300'
            }`}
          />
          {field.state.meta.isTouched && field.state.meta.errors.length > 0 && (
            <p className="text-red-500 text-sm">{field.state.meta.errors[0]}</p>
          )}
        </div>
      )}
    </form.Field>
  );
}
```
