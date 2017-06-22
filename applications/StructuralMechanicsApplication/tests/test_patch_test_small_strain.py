from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestSmallStrain(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)        
        
    
    def _apply_BCs(self,mp,A,b):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(3)
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector()
            u = A*xvec
            u += b
            
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,u)

    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)
        
        g = [0,0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)
        
        if(dim == 2):
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        else:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl) 
            
    def _define_movement(self,dim):
        if(dim == 2):
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;  A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;  A[1,1] = 0.7e-10; A[1,2] = 0.0
            A[2,1] = 0.0;  A[2,1] = 0.0; A[2,2] = 0.0
                    
            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10     
            b[2] = 0.0
            
        else:
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;   A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;   A[1,1] = 0.7e-10; A[1,2] = 0.1e-10
            A[2,1] = -0.2e-10;  A[2,1] = 0.0;     A[2,2] = -0.3e-10
                    
            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10     
            b[2] = 0.7e-10
            
        
        
        return A,b
        
    def _solve(self,mp):
        
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        
        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(mp, 
                                                                        scheme, 
                                                                        linear_solver, 
                                                                        builder_and_solver, 
                                                                        compute_reactions, 
                                                                        reform_step_dofs, 
                                                                        calculate_norm_dx,
                                                                        move_mesh_flag)
        

        #strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp, 
                                                                        #scheme, 
                                                                        #linear_solver, 
                                                                        #convergence_criterion, 
                                                                        #builder_and_solver, 
                                                                        #max_iters, 
                                                                        #compute_reactions, 
                                                                        #reform_step_dofs, 
                                                                        #move_mesh_flag)
        strategy.SetEchoLevel(0)
        
        strategy.Check()
        strategy.Solve()
        
    
    def _check_results(self,mp,A,b):
        
        ##check that the results are exact on the nodes
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(len(b))
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector(2)
            u = A*xvec
            u += b            
            
            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            self.assertAlmostEqual(d[0], u[0])
            self.assertAlmostEqual(d[1], u[1])
            self.assertAlmostEqual(d[2], u[2])
            
    def _check_outputs(self,mp,A,dim):
    
        E = mp.GetProperties()[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        NU =mp.GetProperties()[1].GetValue(KratosMultiphysics.POISSON_RATIO)
        
        #given the matrix A, the analytic deformation graident is F+I
        F = A
        for i in range(3):
            F[i,i] += 1.0
        
        #here compute the Cauchy green strain tensor
        Etensor = KratosMultiphysics.Matrix(3,3)

        for i in range(3):
            for j in range(3):
                Etensor[i,j] = 0.0
                
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Etensor[i,j] += A[k,i]*A[k,j]
                    
        for i in range(3):
            Etensor[i,i] -= 1.0
                    
        for i in range(3):
            for j in range(3):
                Etensor[i,j] = 0.5*Etensor[i,j]
        
        if(dim == 2):
            #verify strain
            reference_strain = KratosMultiphysics.Vector(3)
            reference_strain[0] = Etensor[0,0]
            reference_strain[1] = Etensor[1,1]
            reference_strain[2] = 2.0*Etensor[0,1]
        else:
            reference_strain = KratosMultiphysics.Vector(6)
            reference_strain[0] = Etensor[0,0]
            reference_strain[1] = Etensor[1,1]
            reference_strain[2] = Etensor[2,2]
            reference_strain[3] = 2.0*Etensor[0,1]
            reference_strain[4] = 2.0*Etensor[1,2]
            reference_strain[5] = 2.0*Etensor[0,2]
            
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, mp.ProcessInfo)
            for strain in out:
                for i in range(len(reference_strain)):
                    self.assertAlmostEqual(reference_strain[i], strain[i])
                    
        #finally compute stress
        if(dim == 2):
            #here assume plane stress
            c1 = E / (1.00 - NU*NU);
            c2 = c1 * NU;
            c3 = 0.5* E / (1 + NU);
            reference_stress = KratosMultiphysics.Vector(3)
            reference_stress[0] = c1*reference_strain[0] + c2 * (reference_strain[1])	;
            reference_stress[1] = c1*reference_strain[1] + c2 * (reference_strain[0])	;
            reference_stress[2] = c3*reference_strain[2];
        else:
            c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
            c2 = c1 * ( 1 - NU );
            c3 = c1 * NU;
            c4 = c1 * 0.5 * ( 1 - 2 * NU );
            reference_stress = KratosMultiphysics.Vector(6)
            reference_stress[0] = c2*reference_strain[0] + c3 * (reference_strain[1] + reference_strain[2])
            reference_stress[1] = c2*reference_strain[1] + c3 * (reference_strain[0] + reference_strain[2])
            reference_stress[2] = c2*reference_strain[2] + c3 * (reference_strain[0] + reference_strain[1])
            reference_stress[3] = c4*reference_strain[3]
            reference_stress[4] = c4*reference_strain[4]
            reference_stress[5] = c4*reference_strain[5]
            
        for elem in mp.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, mp.ProcessInfo)
            for stress in out:
                for i in range(len(reference_stress)):
                    self.assertAlmostEqual(reference_stress[i], stress[i],2)        

    def _test_SmallDisplacementElement_2D_triangle(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4])
                
        #create Element
        mp.CreateNewElement("SmallDisplacementElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D3N", 4, [4,1,5], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
        #checking consistent mass matrix
        M = KratosMultiphysics.Matrix(0,0)
        mp.Elements[1].CalculateMassMatrix(M,mp.ProcessInfo)
        Area = mp.Elements[1].GetArea()
        for i in range(3):
            for j in range(3):
                for k in range(dim):
                    if(i==j):
                        coeff = Area/6.0
                    else:
                        coeff = Area/12.0
                    self.assertAlmostEqual(M[i*dim+k,j*dim+k],coeff)
                    
    def _test_SmallDisplacementElement_2D_quadrilateral(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)
        
        #create nodes
        mp.CreateNewNode(1,0.00,3.00,0.00)
        mp.CreateNewNode(2,1.00,2.25,0.00)
        mp.CreateNewNode(3,0.75,1.00,0.00)
        mp.CreateNewNode(4,2.25,2.00,0.00)
        mp.CreateNewNode(5,0.00,0.00,0.00)
        mp.CreateNewNode(6,3.00,3.00,0.00)
        mp.CreateNewNode(7,2.00,0.75,0.00)
        mp.CreateNewNode(8,3.00,0.00,0.00)
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,8])
                
        #create Element
        mp.CreateNewElement("SmallDisplacementElement2D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 5, [2,1,5,3], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
    def _test_SmallDisplacementElement_3D_tetra(self): 
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)        
        
        #create nodes
        mp.CreateNewNode(1,0.00,3.00,3.00)
        mp.CreateNewNode(2,1.00,2.25,3.00)
        mp.CreateNewNode(3,0.75,1.00,3.00)
        mp.CreateNewNode(4,2.25,2.00,3.00)
        mp.CreateNewNode(5,0.00,3.00,0.00)
        mp.CreateNewNode(6,0.00,0.00,3.00)
        mp.CreateNewNode(7,3.00,3.00,3.00)
        mp.CreateNewNode(8,2.00,0.75,3.00)
        mp.CreateNewNode(9,1.00,2.25,0.00)
        mp.CreateNewNode(10,0.75,1.00,0.00)
        mp.CreateNewNode(11,2.25,2.00,0.00)
        mp.CreateNewNode(12,0.00,0.00,0.00)
        mp.CreateNewNode(13,3.00,3.00,0.00)
        mp.CreateNewNode(14,3.00,0.00,3.00)
        mp.CreateNewNode(15,2.00,0.75,0.00)
        mp.CreateNewNode(16,3.00,0.00,0.00)   
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,5,6,7,12,13,14,16])

        #create Element
        mp.CreateNewElement("SmallDisplacementElement3D4N", 1,[12,10,3,15], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 2,[15,8,3,12], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 3,[6,3,8,12], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 4,[12,6,14,8], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 5,[12,16,15,8], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 6,[16,14,8,12], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 7,[13,7,4,16], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 8,[16,13,11,4], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 9,[14,8,4,16], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 10,[16,15,8,4], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 11,[16,14,7,4], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 12,[15,11,4,16], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 13,[5,1,2,13], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 14,[13,5,9,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 15,[7,4,2,13], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 16,[13,11,4,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 17,[13,7,1,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 18,[11,9,2,13], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 19,[11,4,2,15], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 20,[15,11,9,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 21,[8,3,2,15], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 22,[15,10,3,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 23,[15,8,4,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 24,[10,9,2,15], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 25,[10,3,2,12], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 26,[12,10,9,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 27,[6,1,2,12], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 28,[12,5,1,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 29,[12,6,3,2], mp.GetProperties()[1]) 
        mp.CreateNewElement("SmallDisplacementElement3D4N", 30,[5,9,2,12], mp.GetProperties()[1]) 
        
        A,b = self._define_movement(dim)
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
        
    def _test_SmallDisplacementElement_3D_hexa(self): 
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)        
        
        #create nodes
        mp.CreateNewNode(1, 0.00000,  1.00000,  1.00000)
        mp.CreateNewNode(2, 0.16500,  0.74500,  0.70200)
        mp.CreateNewNode(3, 0.27300,  0.75000,  0.23000)
        mp.CreateNewNode(4, 0.78800,  0.69300,  0.64400)
        mp.CreateNewNode(5, 0.32000,  0.18600,  0.64300)
        mp.CreateNewNode(6, 0.00000,  1.00000,  0.00000)
        mp.CreateNewNode(7, 0.00000,  0.00000,  1.00000)
        mp.CreateNewNode(8, 1.00000,  1.00000,  1.00000)
        mp.CreateNewNode(9, 0.67700,  0.30500,  0.68300)
        mp.CreateNewNode(10, 0.24900,  0.34200,  0.19200)
        mp.CreateNewNode(11, 0.85000,  0.64900,  0.26300)
        mp.CreateNewNode(12, 0.82600,  0.28800,  0.28800)
        mp.CreateNewNode(13, 0.00000,  0.00000,  0.00000)
        mp.CreateNewNode(14, 1.00000,  1.00000,  0.00000)
        mp.CreateNewNode(15, 1.00000,  0.00000,  1.00000)
        mp.CreateNewNode(16, 1.00000,  0.00000,  0.00000)       
        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,6,7,8,13,14,15,16])

        #create Element
        mp.CreateNewElement("SmallDisplacementElement3D8N", 1,[10,5,2,3,13,7,1,6], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 2,[12,9,5,10,16,15,7,13], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 3,[12,11,3,10,9,4,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 4,[9,4,2,5,15,8,1,7], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 5,[4,11,3,2,8,14,6,1], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 6,[11,4,9,12,14,8,15,16], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement3D8N", 7,[11,12,10,3,14,16,13,6], mp.GetProperties()[1])
        
        A,b = self._define_movement(dim)
        self._apply_BCs(bcs,A,b)
        self._solve(mp)
        self._check_results(mp,A,b)
        self._check_outputs(mp,A,dim)
    
    def test_execution(self):
        self._test_SmallDisplacementElement_2D_triangle()
        self._test_SmallDisplacementElement_2D_quadrilateral()
        #self._test_SmallDisplacementElement_3D_tetra()
        self._test_SmallDisplacementElement_3D_hexa()
        
if __name__ == '__main__':
    KratosUnittest.main()