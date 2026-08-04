# Authorship and AI-use disclosure

## Maintainer

**Jizhou Guo**  
ORCID: https://orcid.org/0009-0001-0699-9164  
Homepage: https://aster2024.github.io/

Jizhou Guo initiated the project, selected the problem, supplied prompts and source material, steered iterative runs, organized the generated artifacts, requested further verification, and released the repository for public scrutiny.

He does not claim to have independently derived every mathematical argument or personally written every line of verification code.

## Heavy reliance on GPT-5.6-Sol

The development of these proof candidates relied heavily on **OpenAI GPT-5.6-Sol**. Its contribution went substantially beyond editing or routine code completion. Under repeated prompting and feedback, it generated substantial portions of:

- the proof architecture and mathematical exposition;
- new analytic lemmas and proof attempts;
- exhaustive-search and interval-verification code;
- computational certificates and revisions.

Later cross-checks used **Claude** and **OpenAI GPT-5.6 Thinking** to rerun certificates, write differently organized scans, reproduce high-precision KKT solutions, and inspect constants and implementation details.

## Status limitation

These checks remain AI-assisted. Agreement among models or implementations is useful evidence but is not independent human peer review. The repository therefore describes the results as **computer-assisted proof candidates** and uses the statement **“no concrete error has been found so far”**, not “the problem has been definitively solved.”
