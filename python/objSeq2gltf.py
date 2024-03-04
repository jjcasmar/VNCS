from pygltflib import *
import meshio
import base64
import glob
import numpy as np

dir = "/home/jjcasmar/projects/VNCS_Scenes/2DTutorials/Zigzag/Simulation/"
path = "/home/jjcasmar/projects/VNCS_Scenes/2DTutorials/Zigzag/Simulation/coarse_00000.vtu"
path = "/home/jjcasmar/projects/VNCS_Scenes/2DTutorials/Zigzag/Simulation/fine_00000.vtu"

mesh = meshio.read(path)

points = mesh.points.astype('float32')
triangles = mesh.cells[0][1].astype('uint32')

points_buffer_bytes = points.tobytes()
nPoints = len(points)
nTriangles = len(triangles)
triangles_buffer_bytes = triangles.tobytes()
points_buffer_b64 = base64.b64encode(points_buffer_bytes)
triangles_buffer_b64 = base64.b64encode(triangles_buffer_bytes)

gltf = GLTF2()
scene = Scene()

mesh = Mesh()
primitive = Primitive()
primitive.indices = 1
node = Node()
buffer = Buffer()
bufferView1 = BufferView()
bufferView2 = BufferView()
accessor1 = Accessor()
accessor2 = Accessor()


bufferFile = open("0.bin", 'wb')

# add data
buffer.uri = "0.bin"
bufferFile.write(points_buffer_bytes)
bufferFile.write(triangles_buffer_bytes)
buffer.byteLength = len(points_buffer_bytes) + len(triangles_buffer_bytes)

bufferView1.buffer = 0
bufferView1.byteOffset = 0
bufferView1.byteLength = len(points_buffer_bytes)
bufferView1.target = ARRAY_BUFFER

bufferView2.buffer = 0
bufferView2.byteOffset = len(points_buffer_bytes)
bufferView2.byteLength = len(triangles_buffer_bytes)
bufferView2.target = ELEMENT_ARRAY_BUFFER

accessor1.bufferView = 0
accessor1.byteOffset = 0
accessor1.componentType = FLOAT
accessor1.count = nPoints
accessor1.type = VEC3

accessor2.bufferView = 1
accessor2.byteOffset = 0
accessor2.componentType = UNSIGNED_INT
accessor2.count = 3 * nTriangles
accessor2.type = SCALAR

primitive.attributes.POSITION = 0
node.mesh = 0
scene.nodes = [0]

# assemble into a gltf structure
gltf.scenes.append(scene)
gltf.nodes.append(node)
gltf.bufferViews.append(bufferView1)
gltf.bufferViews.append(bufferView2)
gltf.accessors.append(accessor1)
gltf.accessors.append(accessor2)

gltf.buffers.append(buffer)
gltf.meshes.append(mesh)
gltf.meshes[0].primitives.append(primitive)

gltf.convert_buffers(BufferFormat.BINFILE, True)

bufferFile.close()

# save to file
gltf.save("triangle.gltf")
