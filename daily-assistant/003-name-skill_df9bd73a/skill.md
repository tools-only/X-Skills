---
name: nba-data
description: |
  NBA data via ESPN public endpoints — scores, standings, rosters, schedules, game summaries, play-by-play, win probability, injuries, transactions, futures, depth charts, team/player stats, leaders, and news. Zero config, no API keys.

  Use when: user asks about NBA scores, standings, team rosters, schedules, game stats, box scores, play-by-play, injuries, transactions, betting futures, depth charts, team/player statistics, or NBA news.
  Don't use when: user asks about WNBA (use wnba-data), college basketball (use cbb-data), or other sports.
license: MIT
metadata:
  author: machina-sports
  version: "0.1.0"
---

# NBA Data

## Setup

Before first use, check if the CLI is available:
```bash
which sports-skills || pip install sports-skills
```
If `pip install` fails with a Python version error, the package requires Python 3.10+. Find a compatible Python:
```bash
python3 --version  # check version
# If < 3.10, try: python3.12 -m pip install sports-skills
# On macOS with Homebrew: /opt/homebrew/bin/python3.12 -m pip install sports-skills
```
No API keys required.

## Quick Start

Prefer the CLI — it avoids Python import path issues:
```bash
sports-skills nba get_scoreboard
sports-skills nba get_standings --season=2025
sports-skills nba get_teams
```

## Choosing the Season

Derive the current year from the system prompt's date (e.g., `currentDate: 2026-02-18` → current year is 2026).

- **If the user specifies a season**, use it as-is.
- **If the user says "current", "this season", or doesn't specify**: The NBA season runs October–June. If the current month is October–December, the active season year matches the current year. If January–June, the active season started the previous calendar year (use that year as the season).
- **Never hardcode a season.** Always derive it from the system date.

## Commands

### get_scoreboard
Get live/recent NBA scores.
- `date` (str, optional): Date in YYYY-MM-DD format. Defaults to today.

Returns `events[]` with game info, scores, status, and competitors.

### get_standings
Get NBA standings by conference.
- `season` (int, optional): Season year

Returns `groups[]` with Eastern/Western conferences and team standings including W-L, PCT, GB, streak, home/away/conference records, and PPG.

### get_teams
Get all 30 NBA teams. No parameters.

Returns `teams[]` with id, name, abbreviation, logo, and location.

### get_team_roster
Get full roster for a team.
- `team_id` (str, required): ESPN team ID (e.g., "13" for Lakers)

Returns `athletes[]` with name, position, jersey number, height, weight, experience.

### get_team_schedule
Get schedule for a specific team.
- `team_id` (str, required): ESPN team ID
- `season` (int, optional): Season year

Returns `events[]` with opponent, date, score (if played), and venue.

### get_game_summary
Get detailed box score and scoring plays.
- `event_id` (str, required): ESPN event ID

Returns `game_info`, `competitors`, `boxscore` (stats per player), `scoring_plays`, and `leaders`.

### get_leaders
Get NBA statistical leaders (points, rebounds, assists, etc.).
- `season` (int, optional): Season year

Returns `categories[]` with leader rankings per stat category.

### get_news
Get NBA news articles.
- `team_id` (str, optional): Filter by team

Returns `articles[]` with headline, description, published date, and link.

### get_play_by_play
Get full play-by-play data for a game.
- `event_id` (str, required): ESPN event ID

Returns play-by-play detail including period, clock, team, play description, and scoring plays.

### get_win_probability
Get win probability chart data for a game.
- `event_id` (str, required): ESPN event ID

Returns timestamped home/away win probability percentages throughout the game.

### get_schedule
Get NBA schedule for a specific date or season.
- `date` (str, optional): Date in YYYY-MM-DD format
- `season` (int, optional): Season year (used only if no date provided)

Returns `events[]` for the specified date.

### get_injuries
Get current NBA injury reports across all teams. No parameters.

Returns `teams[]` with per-team injury lists including player name, position, status (Out/Doubtful/Questionable/Day-To-Day), injury type, and detail.

### get_transactions
Get recent NBA transactions (trades, signings, waivers).
- `limit` (int, optional): Max transactions to return. Defaults to 50.

Returns `transactions[]` with date, team, and description.

### get_futures
Get NBA futures/odds markets (Championship winner, MVP, etc.).
- `limit` (int, optional): Max entries per market. Defaults to 25.
- `season_year` (int, optional): Season year. Defaults to current.

Returns `futures[]` with market name and entries (team/player name + odds value).

### get_depth_chart
Get depth chart for a specific team.
- `team_id` (str, required): ESPN team ID

Returns `charts[]` with positional depth and player depth order.

### get_team_stats
Get full team statistical profile for a season.
- `team_id` (str, required): ESPN team ID
- `season_year` (int, optional): Season year. Defaults to current.
- `season_type` (int, optional): 1=preseason, 2=regular (default), 3=postseason.

Returns `categories[]` with detailed stats including value, rank, and per-game averages.

### get_player_stats
Get full player statistical profile for a season.
- `player_id` (str, required): ESPN athlete ID
- `season_year` (int, optional): Season year. Defaults to current.
- `season_type` (int, optional): 1=preseason, 2=regular (default), 3=postseason.

Returns `categories[]` with detailed stats including value, rank, and per-game averages.

## Team IDs (Common)

| Team | ID | Team | ID |
|------|-----|------|-----|
| Atlanta Hawks | 1 | Memphis Grizzlies | 29 |
| Boston Celtics | 2 | Miami Heat | 14 |
| Brooklyn Nets | 17 | Milwaukee Bucks | 15 |
| Charlotte Hornets | 30 | Minnesota Timberwolves | 16 |
| Chicago Bulls | 4 | New Orleans Pelicans | 3 |
| Cleveland Cavaliers | 5 | New York Knicks | 18 |
| Dallas Mavericks | 6 | Oklahoma City Thunder | 25 |
| Denver Nuggets | 7 | Orlando Magic | 19 |
| Detroit Pistons | 8 | Philadelphia 76ers | 20 |
| Golden State Warriors | 9 | Phoenix Suns | 21 |
| Houston Rockets | 10 | Portland Trail Blazers | 22 |
| Indiana Pacers | 11 | Sacramento Kings | 23 |
| LA Clippers | 12 | San Antonio Spurs | 24 |
| Los Angeles Lakers | 13 | Toronto Raptors | 28 |
| Utah Jazz | 26 | Washington Wizards | 27 |

## Examples

**User: "What are today's NBA scores?"**
```bash
sports-skills nba get_scoreboard
```

**User: "Show me the Western Conference standings"**
```bash
sports-skills nba get_standings --season=2025
```
Then filter results for Western Conference.

**User: "Who's on the Lakers roster?"**
```bash
sports-skills nba get_team_roster --team_id=13
```

**User: "Show me the full box score for last night's Celtics game"**
1. Find the event_id from `get_scoreboard --date=YYYY-MM-DD`
2. Call `get_game_summary --event_id=<id>` for full box score

**User: "Who's injured on the Lakers?"**
```bash
sports-skills nba get_injuries
```
Then filter results for Los Angeles Lakers (team_id=13).

**User: "What are the NBA championship odds?"**
```bash
sports-skills nba get_futures --limit=10
```

**User: "Show me LeBron's stats this season"**
```bash
sports-skills nba get_player_stats --player_id=1966
```

**User: "How do the Celtics rank in team stats?"**
```bash
sports-skills nba get_team_stats --team_id=2
```

## Error Handling

When a command fails, **do not surface raw errors to the user**. Instead:
1. Catch silently and try alternatives
2. If team name given instead of ID, use `get_teams` to find the ID first
3. Only report failure with a clean message after exhausting alternatives

## Troubleshooting

- **`sports-skills` command not found**: Run `pip install sports-skills`
- **Team not found**: Use `get_teams` to list all teams and find the correct ID
- **No data for future games**: ESPN only returns data for completed or in-progress games
- **Offseason**: `get_scoreboard` returns 0 events — expected. Use `get_standings` or `get_news` instead.
