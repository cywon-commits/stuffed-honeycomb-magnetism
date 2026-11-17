# External archival data and 27-spin checkpoints

The repository retains 1,046 research source/data text files imported from the HPC core transfer bundle (2026-11-17), together with per-stage code, numerical CSV/JSON, plots in SVG, documentation, and a file inventory. Most binary eigenvectors and historical archive ZIP files are **not** checked into ordinary git.

## Archive and checkpoint sources

The following files were previously supplied for this project. Access may require the owner's Google Drive permissions.

- **18-spin / 27-spin finite-size results and source:** https://drive.google.com/file/d/1-jKXpnnyLzUmo0OyjgFlfCfoxEM_uIye/view
- **27-spin stages 1–2 checkpoints:** https://drive.google.com/file/d/1_lx7QD0fpxuq0SO-nrv8lLSZJCW0twek/view
- **27-spin triangular isotropic benchmark:** https://drive.google.com/file/d/15V5uOv_KRn-b5fuiEBLrV8OS4i_PE5c3/view
- **27-spin ground-state level crossings:** https://drive.google.com/file/d/1h_GiE26iTnfj2FzQ-y5uBct53eR9OqTA/view
- **Stage 3 dynamics with verified poles:** https://drive.google.com/file/d/1C7LNsvKVpVF8XdrXpTBdwqcx2JWoLJ7I/view
- **Stage 3 additional spectra / code / numerical data:** https://drive.google.com/file/d/1jVVNWhu0ZH-AbAongotvKcyAziIjvSww/view
- **Stage 3 checkpoint inputs:** https://drive.google.com/file/d/1fdbYZYrtAQIjdA5_7prC4hYshZB9Px4x/view
- **New eigenvectors at r = 0.2375:** https://drive.google.com/file/d/14TraRUJr-xpCZ7cuo-hzsOi4mJwKI170/view
- **New eigenvectors at r = 0.2:** https://drive.google.com/file/d/1M8RCAbANTtFonmmVu9w_gKUzwRC44AC6/view
- **ABC triangle dimer/chirality ED spectra and data:** https://drive.google.com/file/d/11LOo1s2t6YrR3I0NIX6EW2gv58zXJJSp/view

## Provenance and checksum audit

See `data/SOURCE_INVENTORY.csv` (source path, GitHub path where available, size and SHA-256) for exact scope of committed vs externally stored files.

The original `Stuffed_Honeycomb_HPC_Full_2026-11-17.zip` is an approximately 506 MB assembled archive from this project. It must be retained in durable storage separately; this README does **not** assert that this particular archive has been uploaded to GitHub or Drive.

The repository does not contain the full set of binary 27-spin eigensystem checkpoints, and complete dynamical reruns require downloading the appropriate input archives above.

## Technical cautions

The 36-spin low-|Sz| eigensolver is still to be implemented, and existing high-|Sz| 36-spin pilots must not be mistaken for converged ground states. Finite-size crossing locations are geometry-dependent; the physical triangular K point and 3-site cell momentum coordinates are distinct.
