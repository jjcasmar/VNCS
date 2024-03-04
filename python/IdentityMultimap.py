import Sofa
import SofaRuntime
import SofaTypes

def createScene(node):

    node.gravity = [0,0,0]

    node.addObject("RequiredPlugin", name="/home/jjcasmar/projects/CPFSofaPlugin/build/Debug/lib/libCPFSofaPlugin.so")
    node.addObject("DefaultPipeline")
    node.addObject("BruteForceDetection")
    node.addObject("DefaultContactManager")
    node.addObject("DiscreteIntersection")


    odeSolver = node.addObject("EulerImplicitSolver")
    odeSolver.rayleighStiffness = 0.1
    odeSolver.rayleighMass = 0.1

    cgSolver = node.addObject("CGLinearSolver")
    cgSolver.iterations = 25
    cgSolver.tolerance = 1e-5
    cgSolver.threshold = 1e-5

    object1 = node.addChild("object1")
    object1.addObject("MeshObjLoader", name="loader", filename="base.obj")
    object1.addObject("MechanicalObject", template="Vec3d", positions="@loader.positions", name="dof1")

    object2 = node.addChild("object2")
    object2.addObject("GridMeshCreator", name="loader", resolution=[2, 2], translation=[2, 0, 0])
    object2.addObject("Mesh", src="@loader")
    object2.addObject("MechanicalObject", template="Vec3d", src="@loader", name="dof2")

    concat = node.addChild("concatenation")
    concat.addObject("MechanicalObject", template="Vec3d", name="dofAll", showObject=1)
    concat.addObject("MatrixMultiMapping", template=r"Vec3d,Vec3d", input="@../object1/dof1 @../object2/dof2", output="@dofAll")
    concat.addObject("SphereCollisionModel", radius=0.3, selfCollision=1)
    concat.addObject("UniformMass", vertexMass=1)
    concat.addObject("ConstantForceField", indices=[0], force=[1, 0, 0])
