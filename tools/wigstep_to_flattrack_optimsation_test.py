
"""

This is for optimisng flat_track generation which is really too slow at the moment.

"""
import sys, os
sys.path.append(os.path.realpath("../../"))
import glbase as gl

if __name__ == "__main__":
    # If you want to try this, you'll need to make a wigfix file
    # Grab about ~4 Mbp of wigfix file, preferably from a few different chromosomes.
    # Optimisation testing.
    import cProfile
    import pstats
    
    cProfile.run("gl.wigstep_to_flat('test_wigfix.wig', 'test_flat.flat', name='meh', bin_format='f')", "profile_out")
    p = pstats.Stats("profile_out")
    p.strip_dirs().sort_stats('cumulative').print_stats()
    
    flat = gl.flat_track(filename='test_flat.flat', name='meh', bin_format='f')
    print(flat.get("chr1:3000306-3000326"))
    
    # Correct answer is:
    """
    array('f', [0.22300000488758087, 0.23999999463558197, 0.22300000488758087, 0.26199999451637268, 
    0.22300000488758087, 0.23600000143051147, 0.26199999451637268, 0.22300000488758087, 
    0.22300000488758087, -1.3509999513626099, 0.23999999463558197, 0.22300000488758087, 
    -0.9089999794960022, 0.22300000488758087, 0.22300000488758087, 0.22300000488758087, 
    0.26199999451637268, 0.23999999463558197, 0.22300000488758087, 0.22300000488758087, 
    0.23600000143051147])
    """
    
    #print flat.get("chr1:3159900-3160100")