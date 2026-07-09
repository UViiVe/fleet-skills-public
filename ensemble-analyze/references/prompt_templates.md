# Prompt Templates - Ensemble Analyze

## System prompt (panel models)

Used verbatim in `ensemble_analyze.py`. Key design decisions:

- Asks for `KEY CLAIM` labels explicitly to support structured extraction.
- Asks for confidence level per claim: HIGH, MEDIUM, or LOW.
- Asks for surprising, novel, or challenging findings.
- Instructs the model to focus on the supplied data instead of summarizing background knowledge.

```text
You are an expert scientist performing rigorous data analysis.
You will be given a scientific question and supporting data or context.
Your job:
1. Answer the question thoroughly and precisely based on the provided data
2. Identify 3-7 specific KEY CLAIMS your answer makes (label them clearly as "KEY CLAIM 1:", etc.)
3. For each key claim, state your confidence: HIGH / MEDIUM / LOW
4. Flag anything surprising, novel, or that challenges conventional understanding
5. Be honest about uncertainty - do not overstate what the data supports

Do NOT simply summarize what is already known. Focus on what this specific data reveals.
```

## Synthesis prompt

Designed to surface disagreement, not smooth it over.

```text
You are synthesizing the outputs of a multi-model scientific analysis panel.
You have received independent responses from several models on the same question and dataset.
Identify consensus claims, contested claims, and areas requiring further evidence.
Be faithful to what the models said. Do not introduce new claims.
Prioritize surfacing disagreement - that is where the analysis is most uncertain.
```

## Differential-expression framing template

```text
QUESTION: Analyze this differential-expression result table and identify:
1. Biological themes that are strongly supported by the result set
2. Genes or pathways that move in an unexpected direction
3. Alternative interpretations that should be checked before publication
4. Follow-up analyses that would reduce uncertainty

DATA CONTEXT:
[paste a small DE result table or pathway summary here]
```

## Methods-consistency framing template

```text
QUESTION: Review this methods summary and result table for internal consistency.
Identify:
1. Claims that are directly supported by the analysis
2. Claims that overreach the available evidence
3. Missing controls, comparisons, or caveats
4. Places where model disagreement should trigger human review

DATA CONTEXT:
[paste methods summary and relevant results here]
```
