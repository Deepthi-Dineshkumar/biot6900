# BIOT 6900 · Module 4: Structure-Based Druggability of Your Targets

**Name:** Deepthi Dineshkumar

**Disease:** Pancreatic ductal adenocarcinoma (PDAC) — CPTAC-PDAC cohort via LinkedOmics, carried forward from Modules 2 and 3

**Targets:**

| Gene | UniProt | Module 2 rank | PDB used | Call |
|---|---|---|---|---|
| CEL — carboxyl ester lipase | P19835 | 4 | 1F6W chain A | Possible but uncertain |
| XDH — xanthine dehydrogenase/oxidase | P47989 | 9 | 2E1Q chain A | Structure supports a small-molecule pocket |

## What didn't work

- **fpocket would not install on Windows.**Step 6 was run on the Northeastern Explorer cluster: the two AlphaFold models   were uploaded 
  `fpocket -f` was run on each, and the `AF-P19835_out` and `AF-P47989_out` folders were copied
  back into `structures/`. The Step 6 cell was modified to reuse those folders if they exist, so
  the notebook still runs top to bottom on Windows.
- **No ligand-site check for CEL.** 1F6W contains no eligible bound small molecule, so the
  ligand-site check printed a one-line message for CEL and ran only for XDH. This is the main
  reason CEL's call is "possible but uncertain" rather than "supports a small-molecule pocket".
- **No inter-domain PAE for either target.** The domain finder returned a single confident group
  per protein (CEL 21-552, XDH 2-1333), so `inter_pae` is NaN in both rows of the summary CSV.
  That is the expected output for a single group, not a failure.
- **2E1Q is the E803V mutant** of human XOR. Because the Step 5 superposition keeps only pairs
  where the amino acid is identical in both structures, position 803 is excluded automatically.
- **Numbering offset.** 1F6W uses mature-protein numbering, offset by the 20-residue signal
  peptide from UniProt numbering. All positions reported are UniProt numbering; the
  Needleman-Wunsch alignment in Step 5 handles the conversion.
