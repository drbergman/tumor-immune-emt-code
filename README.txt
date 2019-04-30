We develop an agent-based model to describe the relationships between cancer, the immune system, and EMT, building on cell-cycle and tissue-cell components of a previous model. 
The agents in the model are the tissue cells with tumorigenic potential.
These can either be mutation-free or, resulting from DNA damage during the cell cycle, can have any combination of three possible pathway mutations.
A single mutation is sufficient for the cell to be considered as a mutant cell.

We model immune cells as continuous variables, appropriate since at the level of immune cells it is not necessary to keep track of single cell effects.
The cytokine TGFB is continuous within the TME, i.e. the system is well-mixed.
Tissue cells can take on either epithelial or mesenchymal phenotypes in a plastic manner: these phenotypes depend on both the TME and cell-intrinsic factors.
While the score is continuous, a threshold determines if a given cell is labeled as epithelial or mesenchymal.
Even though each cell does have an EMT score on a spectrum of values, only the labeling influences its interactions in the model.
