# Daily Reddit Workflow for fabrikIQ

**Account:** u/Ok-Painter2695
**Frequenz:** Täglich, 10-15 Minuten

---

## Quick-Start Befehl

```
Scanne Reddit nach Posts für fabrikIQ. Zeige relevante Posts, bereite Kommentare vor, warte auf meine Freigabe vor dem Posten.
```

---

## Workflow-Schritte

### 1. Reddit Scan (automatisch)

**WICHTIG: Nur Bash + curl verwenden!**
- WebFetch, Ref, Exa → blockiert von Reddit
- Playwright/Chrome → funktioniert, aber langsam
- **curl mit Browser User-Agent → funktioniert zuverlässig**

```bash
# Diese Subreddits werden gescannt (Browser User-Agent PFLICHT!):
curl -s -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36" "https://www.reddit.com/r/manufacturing/new.json?limit=15"
curl -s -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36" "https://www.reddit.com/r/PLC/new.json?limit=10"
curl -s -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36" "https://www.reddit.com/r/LeanManufacturing/new.json?limit=10"
```

**Parsing mit Python:**
```bash
curl -s -A "Mozilla/5.0..." "https://www.reddit.com/r/manufacturing/new.json?limit=10" | python -c "
import json, sys
data = json.load(sys.stdin)
for post in data['data']['children']:
    p = post['data']
    print(f\"Title: {p['title'][:60]}\")
    print(f\"Author: u/{p['author']} | Score: {p['score']} | Comments: {p['num_comments']}\")
    print(f\"URL: https://reddit.com{p['permalink']}\")
    print('---')
"
```

### 2. Filterung nach Keywords
- MES, OEE, SPC, Quality, Measurement
- Data analysis, Excel problems
- Small manufacturer, automation
- AI in manufacturing

### 3. Relevanz-Bewertung
| Score | Kriterien |
|-------|-----------|
| HOCH | Direkt fabrikIQ-Thema (Datenanalyse, OEE, Qualität) |
| MITTEL | Manufacturing allgemein, könnte helfen |
| NIEDRIG | Off-topic oder nur News |

### 4. Kommentar-Vorbereitung
**Stil-Regeln:**
- Alles kleingeschrieben bei kurzen Kommentaren
- Variierende Satzlängen (kurz. dann länger erklärend.)
- Kleine Fehler einbauen (wasnt, theres, teh, adn, wierd)
- Emotional, authentisch, nicht perfekt
- Erwähne Zwillinge nur wenn es passt
- NIEMALS direkte Produktwerbung
- KEINE führenden Leerzeichen vor Sätzen (typisch KI!)
- KEINE Markdown-Formatierung (keine **Sterne** für Bold, keine Listen mit -)
- Überschriften als Plaintext: "Gatekeeper limits -" statt "**Gatekeeper limits** -"
- Abrupte Enden OK - nicht jeder Kommentar braucht Fazit
- Reddit-Slang: tbh, imo, ngl, idk, gonna, kinda
- Doppeltes Leerzeichen gelegentlich OK

### 5. Freigabe-Loop
1. Zeige Post-Titel + URL
2. Zeige vorgeschlagenen Kommentar
3. Warte auf "ja" oder Änderungswünsche
4. User postet selbst (Browser)
5. Weiter zum nächsten

---

## Tracking in Marketing.md

Nach jeder Session aktualisieren:
- Welche Posts kommentiert
- Welche Insights für ROADMAP
- Welche User für Beta-Tester Outreach

---

## Beispiel-Session

```
User: "Scanne Reddit für fabrikIQ"

Claude:
## Gefundene Posts (Top 3)
1. "Help with OEE tracking" - r/manufacturing - 12 comments - HOCH
2. "Excel vs MES debate" - r/manufacturing - 8 comments - HOCH
3. "PLC data logging" - r/PLC - 3 comments - MITTEL

## Kommentar 1
[Post-Details]
[Vorgeschlagener Text]

Freigabe? (ja/ändern)

User: "ja"

Claude: Gepostet! Weiter zu Kommentar 2...
```

---

## Zeitplan-Empfehlung

| Tag | Zeit (CET) | Fokus |
|-----|------------|-------|
| Mo | 9:00 | Scan + 2-3 Kommentare |
| Di | 17:00 | Check Antworten, Follow-up |
| Mi | 9:00 | Scan + 2-3 Kommentare |
| Do | 17:00 | Eigener Post vorbereiten |
| Fr | 9:00 | Leichter Scan, Wochenend-Pause |

**Beste Zeiten für USA-Reichweite:** 15:00-17:00 CET (9-11 AM EST)

---

## Metriken tracken

- Kommentare pro Woche: Ziel 10+
- Upvotes auf Kommentare
- DMs/Anfragen von interessierten Usern
- Karma-Entwicklung auf u/Ok-Painter2695
