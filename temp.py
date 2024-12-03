from chimerax_utils import * 
c = ChimeraCommandManager()
gt_path = "/Users/nadavsellam/Documents/research/proteinx_guidance/pdbs/4ole_A.pdb"
generated_path = "/Users/nadavsellam/Documents/research/proteinx_guidance/output/4ole/seed_101/predictions/4ole_sample_0.cif"

def load_scene():
    c.clear()
    gt = c.load_protein(gt_path)
    generated = c.load_protein(generated_path)
    generated.align(gt)
load_scene()