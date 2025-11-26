from slophep.Generator import MCGenerator

try:
    import ROOT
    import array

    def generate_to_root(gen: MCGenerator, 
                         N: int, 
                         outpath: str, 
                         outtree: str = "DecayTree",
                         printprog: int = 100):
        if gen.maxBF < 0:
            print("Probing for max BF for reject method!")
            gen.get_max_BF()
            print(f"Obtained maxBF = {gen.maxBF}")

        coordarr = {}
        outfile = ROOT.TFile(outpath, "RECREATE")
        tree = ROOT.TTree(outtree, outtree)
        for icoord in gen.coords:
            coordarr[icoord] = array.array('d', [0])
            tree.Branch(icoord, coordarr[icoord], f'{icoord}/D')
        
        for i in range(N):
            if i % abs(printprog) == 0:
                print(f"Generating event {i}")
            point = gen.generate_point()
            for iname, ival in zip(gen.coords, point):
                coordarr[iname][0] = ival
            tree.Fill()
        
        outfile.Write()
        outfile.Close()
        return

except:
    print("ROOT could not be imported, cannot use interface for MCGenerator!")