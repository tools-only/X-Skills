# Book Metrics Skill Creation

Use the skill-creator skill to create a new skill called book-metrics-generator.
This skill will execute a single shell script in @scripts/book-metrics-generator.sh

This shell script will create two files.  

1. The first file is called called @docs/learning-graph/book-metrics.md. It contains 
all the overall metrics about the book including things like chapter count, word 
count, faq count, concept count etc.  The one thing that the book-metrics.md file does
 NOT do is store metrics about individual chapters.

2. The second file that the @script/book-metrics-generator.sh does is create a second 
markdown file called @docs/learning-graph/chapter-metric.md.  This file only contains 
chapter-by-chapter metrics, one row in a Markdown table per chapter.

The UNIX shell script is just a small wrapper that calls a powerful python program 
that does all the work.

The metrics files will actually be created by a python program you create in 
@src/book-metrics/book-metrics.py

The book-metrics.md file will contain a four column markdown table with the following 
collumns:  1. Metric Name 2. Metric Value 3. Link (when appropriate). 4. Notes

Here are the following book metrics:

1. Chapters
2. Concepts
3. Glossary Terms
4. FAQs
5. Quiz Questions
6. Diagrams (level 4 headers that start with "#### Diagram:"
7. Equations (LaTeX expressions)
8. MicroSims (count of the directories in the @doc/sims folder)
9. Total Words
10. Equivalent Pages (make resonable assumptions on words per page and the space taken
 by diagrams)

If there are any other metrics you think are releven (like number of links) feel free 
to include these.

The table will be created by a python program called from 
@scripts/book-metrics-generator.sh that you will create.

The chapter-metrics.md file is also created by the python program.  It creates a 
Markdown table with the following columns:

1. Chapter Number (with leading zeros removed)
2. Chapter Name
3. Sections (count of markdown level 2 and 3 headers)
4. Diagrams
5. Word Count

Text after the table explains what is being counted.

The python program is modular and flexible and it is easy to add new book or chapter 
level metrics. 

⏺ I'll help you create the book-metrics-generator skill using the skill-creator skill.
   This will generate comprehensive metrics about intelligent textbooks.

