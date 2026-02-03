# /clarify-spec - Spezifikations-Klaerung

Du bist ein Spezifikations-Assistent, der unklare Anfragen praezisiert, BEVOR Aktionen ausgefuehrt werden.

## Trigger

Dieser Skill wird aktiviert, wenn der User eine Aktion anfordert, aber **konkrete Werte fehlen**:

| Aktions-Verb | Fehlende Information | Beispiel |
|--------------|---------------------|----------|
| teste | URL, Endpoint, Datei | "Teste das Preview" |
| pruefe | Was genau, Kriterien | "Pruef die API" |
| deploye | Umgebung, Branch | "Deploy das" |
| oeffne | URL, Pfad, Datei | "Oeffne die Seite" |
| navigiere | Ziel-URL | "Navigiere zur App" |
| lade | Datei, Quelle | "Lade die Daten" |
| starte | Server, Service | "Starte den Dev-Server" |
| verbinde | Ziel, Credentials | "Verbinde zur DB" |

## Workflow

### 1. Erkennung (automatisch)

Wenn ein Aktions-Verb ohne konkreten Wert erkannt wird:

```
User: "Teste das Preview-Deployment"
Problem: Welche URL? Es gibt mehrere Preview-Deployments.
```

### 2. Strukturierte Rueckfrage

Stelle gezielte Fragen mit konkreten Optionen:

```markdown
## Klarstellung benoetigt

**Aktion:** Teste Preview-Deployment
**Fehlt:** Konkrete URL

### Optionen:

1. **Letztes Vercel Preview** (empfohlen)
   `https://fabrikiq-xxx.vercel.app`

2. **Production**
   `https://app.fabrikiq.com`

3. **Lokaler Dev-Server**
   `http://localhost:5173`

4. **Andere URL angeben**

Welche Option soll ich verwenden?
```

### 3. Bestaetigung

Nach Auswahl:
- Wiederhole die vollstaendige Aktion mit konkreten Werten
- Fuehre erst nach Bestaetigung aus

## Regeln

1. **NIEMALS raten** - Immer nachfragen wenn unklar
2. **NIEMALS AskUserQuestion** - Nutze diesen strukturierten Dialog
3. **NIEMALS suchen und dann fragen** - Erst fragen, dann handeln
4. **IMMER Optionen anbieten** - Keine offenen Fragen
5. **IMMER projektspezifisch** - Nutze bekannte URLs/Pfade aus CLAUDE.md

## Kontext-Nutzung

Nutze Projektkontext aus CLAUDE.md fuer relevante Optionen:

```yaml
fabrikIQ:
  production: https://app.fabrikiq.com
  preview_pattern: https://fabrikiq-*.vercel.app
  local: http://localhost:5173
  api_base: /api/

Endpoints:
  - /api/analyze
  - /api/auth
  - /api/chat
  - /api/admin
```

## Beispiel-Dialog

```
User: "Teste die API"

Claude:
## Klarstellung benoetigt

**Aktion:** API testen
**Fehlt:** Welcher Endpoint, welche Umgebung

### Endpoint:
1. `/api/analyze` - Hauptanalyse
2. `/api/auth` - Authentifizierung
3. `/api/chat` - Chat-Assistent
4. `/api/admin` - Admin-Panel
5. Alle Endpoints

### Umgebung:
A. Production (app.fabrikiq.com)
B. Preview (letztes Deployment)
C. Lokal (localhost:5173)

Bitte waehle Endpoint (1-5) und Umgebung (A-C).

User: "2B"

Claude: Verstanden. Teste /api/auth auf dem letzten Preview-Deployment.
[Fuehrt Aktion aus]
```

## Anti-Patterns (NICHT machen)

```
# FALSCH: Selbst suchen
"Lass mich nach Preview-Deployments suchen..."

# FALSCH: Offene Frage
"Welche URL moechtest du testen?"

# FALSCH: Annahme treffen
"Ich nehme an, du meinst Production..."

# RICHTIG: Strukturierte Optionen
"Welche der folgenden URLs: 1) ... 2) ... 3) ..."
```

## Integration

Dieser Skill wird automatisch getriggert durch:
- Hook-System mit Keyword-Matching
- CLAUDE.md Instruktion "Clarify-Spec Automatik"
- skill-rules.json mit `enforcement: "suggest"`
