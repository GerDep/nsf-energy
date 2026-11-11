import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import numpy as np

from analysis.reorganize_variance import fold22, variance, summarize
from algebra import bch
from algebra import extended_weights as weights
from scripts.import_reference_data import install, verified_files


class FoldedMoments(unittest.TestCase):
    def test_reciprocal_pair_matches_explicit_energy_distribution(self):
        # Energies 2 and 1/2, weights 1 and 4.
        m=fold22([dict(u=2.,P=5.,D=-1.)])
        np.testing.assert_allclose(m,[5.,-1.,2.,5.],rtol=0,atol=1e-14)
        self.assertAlmostEqual(variance(m),.36)

    def test_upper_and_lower_only_branches_match_explicit_energies(self):
        m=fold22([dict(u=3.,P=2.,D=4.)],maximum=3.)
        np.testing.assert_allclose(m,[2.,4.,8.,18.],rtol=0,atol=1e-14)
        m=fold22([dict(u=10.,P=3.,D=-2.7)],maximum=2.)
        np.testing.assert_allclose(m,[3.,-2.7,2.43,.03],rtol=0,atol=1e-14)

    def test_total_variance_and_covariance_identity(self):
        # Channel 2 has equally weighted E=1/4 and 3/4. The channel-3 replicas
        # contain symmetric two-point distributions, of changing width.
        m2=np.array([1.,-.5,.3125,.3125])
        g=np.array([[2.,1.,1.,5.],[2.,1.2,1.22,5.62],[2.,.8,.82,4.42]])
        s=summarize(m2,[g])
        E=np.array([.25,.75,1.,2.,1.1,2.1,.9,1.9])
        w=np.array([.5,.5]+[1/3]*6)
        mu=np.sum(w*E)/np.sum(w)
        expected=np.sum(w*(E-mu)**2)/np.sum(w)
        self.assertAlmostEqual(s['v'],expected,places=14)
        self.assertAlmostEqual(s['v'],s['within']+s['between'],places=14)
        self.assertLess(s['checks']['raw_vs_centered_variance'],1e-12)
        self.assertGreater(s['v_SE'],0)


class BCHAndCoherence(unittest.TestCase):
    def test_nonzero_quintic_is_reconstructed(self):
        def hermitian(cre,ann,z):return {(cre,ann):z,(ann,cre):z.conjugate()}
        T1=hermitian((2,),(0,1),.3+.2j)
        T2=hermitian((1,2,2),(0,),-.1+.15j)
        T3=hermitian((1,2,2),(0,1),.17+.09j)
        out,_,_=bch.outgoing4(T1,T2,T3,3)
        rebuilt,check=bch.reconstruct_quintic(out,T1,T2,3)
        self.assertLess(bch.norm(bch.add(rebuilt,bch.scale(T3,-1))),1e-12)
        self.assertLess(check['hermiticity_residual'],1e-12)

    def test_missing_physical_coefficient_is_rejected(self):
        with self.assertRaises(ValueError):bch.reconstruct_quintic(None,{}, {},3)
        a=np.zeros((1,8,4),complex)
        with self.assertRaises(weights.MissingPhysicalInput):
            weights.channel23(a,a,a,None,[2.],[1.],source_record={'scope':'synthetic test'})

    def test_interference_is_coherent(self):
        vals=[.3+.2j,.1-.4j,-.15+.05j,-.2+.3j]
        arrays=[]
        for z in vals:
            a=np.zeros((1,8,4),complex);a[0,0,0]=z;arrays.append(a)
        r=weights.channel23(*arrays,[2.],[2.],source_record={'scope':'synthetic algebra test'})
        expected=2*abs(sum(vals))**2/4
        self.assertAlmostEqual(r['orders'][6]['P'],expected,places=14)
        self.assertAlmostEqual(sum(v['P'] for v in r['six_blocks'].values()),expected,places=14)


class DataImport(unittest.TestCase):
    def test_manifest_round_trip_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);payload=b'{"sample":1}\n'
            meta={'schema':'nsf-reference-data-v1','files':{'sample.json':{'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}}}
            path=p/'data.zip'
            with zipfile.ZipFile(path,'w') as z:
                z.writestr('sample.json',payload);z.writestr('manifest.json',json.dumps(meta))
            self.assertEqual(install(path,p/'reference'),1)
            self.assertEqual((p/'reference/sample.json').read_bytes(),payload)
            self.assertEqual(install(path,p/'reference'),1)
            meta['files']['sample.json']['sha256']='0'*64
            with zipfile.ZipFile(path,'w') as z:
                z.writestr('sample.json',payload);z.writestr('manifest.json',json.dumps(meta))
            with self.assertRaises(ValueError):verified_files(path)


class PortableDriver(unittest.TestCase):
    def test_embedded_engine_runs_without_workspace_dependencies(self):
        import graficos_nsf_local as driver
        with tempfile.TemporaryDirectory(prefix='nsf driver ') as tmp:
            root,_=driver.install_engine(Path(tmp))
            driver.init_worker(str(root))
            point=driver._MODEL.evaluate22((1.2,.3,4,1.,1.,0,0))
            self.assertTrue(np.isfinite(point['weight']))
            self.assertGreater(point['weight'],0)


if __name__=='__main__':unittest.main()
