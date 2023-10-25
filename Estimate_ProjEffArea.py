import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mplcursors
import pyvista as pv
import math
from matplotlib.patches import Polygon

def estimate_projeffarea(height_incd, vertices_incd, height_port_incd, center_port, radius_port, height_recv_port):
    """
    输入：
    height_incd: 入射孔径的平面高度
    vertices_incd: 入射孔径顶点在平面内的（二维）位置坐标
    height_port_incd: 接口孔径到入射孔径的平面距离（正值）
    center_port: 接口孔径圆心在平面内的（二维）位置坐标
    radius_port: 接口孔径圆的半径
    height_recv_port: 接收板到接口孔径的平面距离（正值）
    输出：
    接收板二维格点的标记阵列 
    视场二维格点的标记阵列
    视场在全天空的Mollweide投影的标记阵列
    """
    # 入射孔径内切圆圆心的轨迹（加上以其边缘为圆心的外包络，构成视场的full view）
    r_incd_proj = radius_port * (height_port_incd + height_recv_port) / height_recv_port
    inner_vertices = calculate_inner_vertices(vertices_incd, r_incd_proj)
    for inner_vertex in inner_vertices:
        inner_vertex = list(inner_vertex) + [height_incd]
        inner_vertex = tuple(inner_vertex)
    # 入射孔径内切圆圆心的轨迹，以接口孔径圆心为投影点，在接收板上的投影（构成接收板的全影区域）
    projection_point = list(center_port) + [height_incd - height_port_incd]
    projection_point = tuple(projection_point)
    full_shadow_vertices = calculate_projection(inner_vertices, height_incd, projection_point, height_recv_port)
    # 入射孔径以接口孔径为投影点，在接收板上的投影（加上以其边缘为圆心的外包络，构成接收板的半影区域）
    for vertex_incd in vertices_incd:
        vertex_incd = list(vertex_incd) + [height_incd]
        vertex_incd = tuple(vertex_incd)
    part_shadow_vertices = calculate_projection(vertices_incd, height_incd, projection_point,height_recv_port)
    # 绘制接收板二维格点的标记阵列
    r_recv_proj = radius_port * (height_port_incd + height_recv_port) / height_port_incd
    plot_recv_distribution(full_shadow_vertices, part_shadow_vertices, r_recv_proj)
    # 绘制视场二维格点的标记阵列、视场在全天空的Mollweide投影的标记阵列
    plot_incd_distribution(vertices_incd, inner_vertices, r_incd_proj, center_port, height_port_incd)
    plt.show()
    return

def calculate_projection(vertices, vertices_height, projection_point, projected_height):
    """
    输入：
    vertices: 作为投影的“物”的顶点的（二维）位置坐标
    vertices_height: 作为投影的“物”的顶点的高度
    projection_point: 投影点的（三维）位置坐标
    projection_height: 作为投影的“像”的平面，距离投影点的高度差（正值）
    输出：
    projected_vertices: 投影的“像”的顶点的（三维）位置坐标
    """
    vertices = np.array(vertices)
    projection_point = np.array(projection_point)
    projected_vertices = []
    for vertex in vertices:
        vertex = list(vertex) + [vertices_height]
        direction = projection_point - vertex
        projected_vertex = projection_point + direction * projected_height / np.abs(direction[2])
        projected_vertices.append(projected_vertex)
    return projected_vertices

def calculate_inner_vertices(vertices, r):
    """
    输入：
    vertices: 多边形顶点的（二维）坐标
    r: 内切圆的半径
    输出：
    inner_vertices: 内切圆圆心轨迹的（二维）坐标
    """
    inner_vertices = []
    n = len(vertices)
    
    for i in range(n):
        p1 = vertices[i]
        p2 = vertices[(i + 1) % n]
        p3 = vertices[(i + 2) % n]
        # 计算内角平分线的矢量方向并归一化
        edge1 = (p2[0] - p1[0], p2[1] - p1[1])
        edge2 = (p2[0] - p3[0], p2[1] - p3[1])
        length1 = np.sqrt(edge1[0] ** 2 + edge1[1] ** 2)
        length2 = np.sqrt(edge2[0] ** 2 + edge2[1] ** 2)
        edge1 = (edge1[0] / length1, edge1[1] / length1)
        edge2 = (edge2[0] / length2, edge2[1] / length2)
        # 计算内角半角的正弦值
        angle_half_sin = edge1[0] * edge2[1] - edge1[1] * edge2[0]
        # 计算新顶点的位置
        if abs(angle_half_sin) > 1e-6:
            length = r / angle_half_sin
            inner_vertex_x = p2[0] + length * (edge1[0] + edge2[0])
            inner_vertex_y = p2[1] + length * (edge1[1] + edge2[1])
            inner_vertex = (inner_vertex_x, inner_vertex_y)
            inner_vertices.append(inner_vertex)
    return inner_vertices

def point_polygon_relationship(point, polygon, r):
    """
    输入：
    point: 待判断的点的（二维）位置坐标
    polygon: 多边形顶点的（二维）位置坐标
    r: 外包络圆的半径
    输出：
    flag=2: 点位于多边形内部
    flag=1: 点位于多边形外部，且到多边形的距离小于r
    flag=0: 点到多边形的距离总是大于r
    """

    # 判断点是否在多边形内部
    def is_inside_polygon(p, poly):
        n = len(poly)
        inside = False
        x, y = p

        p1x, p1y = poly[0]
        for i in range(1, n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside
    if is_inside_polygon(point, polygon):
        return 2

    # 判断点到多边形顶点的距离
    min_distance_to_vertex = min(math.sqrt((point[0] - v[0]) ** 2 + (point[1] - v[1]) ** 2) for v in polygon)
    if min_distance_to_vertex <= r:
        return 1

    # 判断点到多边形的边的距离
    for p1, p2 in zip(polygon, polygon[1:] + [polygon[0]]):
        x1, y1 = p1
        x2, y2 = p2
        # 计算多边形边所在直线的距离
        distance_to_line = abs((x2 - x1) * (y1 - point[1]) - (x1 - point[0]) * (y2 - y1)) / math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if distance_to_line <= r:
            # 计算点到边两端顶点的矢量夹角
            vector1 = (x1 - point[0], y1 - point[1])
            vector2 = (x2 - point[0], y2 - point[1])
            dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1]
            length1 = math.sqrt(vector1[0] ** 2 + vector1[1] ** 2)
            length2 = math.sqrt(vector2[0] ** 2 + vector2[1] ** 2)
            if dot_product / (length1 * length2) < math.cos(math.atan(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / r)):
                # 如果点到边两端顶点的矢量夹角大于临界情况（正切值为边长/r），则点到边的距离满足条件；取其cos值则为小于
                return 1

    # 点到多边形任何一点的距离均不满足条件
    return 0

def plot_incd_distribution(vertices_incd, inner_vertices, r_incd_proj, center_port, height_port_incd):
    """
    输入：
    vertices_incd: 入射孔径顶点在平面内的（二维）位置坐标
    inner_vertices: 内切圆圆心轨迹的（二维）坐标，加上以其边缘为圆心的外包络，构成视场的full view
    r_incd_proj: 从接收板看过去，接口孔径在入射平面内的半径
    center_port: 接口孔径圆心在平面内的（二维）位置坐标
    height_port_incd: 接口孔径到入射孔径的平面距离（正值）
    输出：
    视场二维格点的标记阵列
    视场在全天空的经纬度坐标的标记阵列
    视场在全天空的Mollweide投影的标记阵列
    """
    # 生成二维网格
    x_min, x_max, y_min, y_max = -2, 6, -2, 6
    grid_size = 0.05
    xx, yy = np.meshgrid(np.arange(x_min, x_max, grid_size), np.arange(y_min, y_max, grid_size))
    lonn, latt = xx, yy

    # 标记网格点位置关系
    marked_points = np.zeros_like(xx, dtype=int)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            point = (xx[i, j], yy[i, j])
            if  point_polygon_relationship(point,vertices_incd,r=0.01) == 0:
                marked_points[i, j] = 0
            elif point_polygon_relationship(point,inner_vertices,r_incd_proj) == 2 \
                or point_polygon_relationship(point,inner_vertices,r_incd_proj) == 1:
                marked_points[i, j] = 2
            else:
                marked_points[i, j] = 1

    # 绘制标记结果
    cmap = plt.cm.get_cmap('coolwarm', 3)
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(xx, yy, marked_points, cmap=cmap, shading='auto')
    cbar = plt.colorbar()
    cbar.set_ticks([0, 1, 2])
    cbar_labels = ['No View', 'Partial View', 'Full View']
    cbar.set_ticklabels(cbar_labels)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.grid(color='k')

    # 绘制多边形
    polygon_patch_full = plt.Polygon(vertices_incd, closed=True, fill=False, edgecolor='blue', lw=2)
    polygon_patch_part = plt.Polygon(inner_vertices, closed=True, fill=False, edgecolor='green', lw=2)
    plt.gca().add_patch(polygon_patch_full)
    plt.gca().add_patch(polygon_patch_part)
    plt.title('Distribution in Incident Plane')
    
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            dx = xx[i ,j] - center_port[0]
            dy = yy[i, j] - center_port[1]
            dz = height_port_incd
            lonn[i, j] = np.arctan2(dx, dz)
            latt[i, j] = np.arctan2(dy, dz)
                
    # 绘制经纬度坐标的标记阵列
    plt.figure(figsize=(8, 4)) 
    plt.pcolormesh(lonn, latt, marked_points, cmap=cmap, shading='auto')
    cbar = plt.colorbar()
    cbar.set_ticks([0, 1, 2])
    cbar_labels = ['No View', 'Partial View', 'Full View']
    cbar.set_ticklabels(cbar_labels)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Distribution in Lon-Lat Plane')
    plt.grid(color='k')
    
    # 绘制Mollweide投影的标记阵列
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot()
    ax.patch.set_facecolor('k')
    x_mollweide = 2 * np.sqrt(2) * np.cos(latt) * np.sin(lonn / 2)
    y_mollweide = np.sqrt(2) * np.sin(latt)
    plt.pcolormesh(x_mollweide, y_mollweide, marked_points, cmap=cmap, shading='auto')
    cbar = plt.colorbar()
    cbar.set_ticks([0, 1, 2])
    cbar_labels = ['No View', 'Partial View', 'Full View']
    cbar.set_ticklabels(cbar_labels)
    plt.xlabel('Mollweide X')
    plt.ylabel('Mollweide Y')
    plt.title('Distribution in Mollweide Projection')
    
    # 在Mollweide投影图中添加经纬线
    lon_line = np.linspace(np.min(lonn), np.max(lonn), 21)
    lat_line = np.linspace(np.min(latt), np.max(latt), 21)
    lon_grid, lat_grid = np.meshgrid(lon_line, lat_line)
    x_mollweide_grid = 2 * np.sqrt(2) * np.cos(lat_grid) * np.sin(lon_grid / 2)
    y_mollweide_grid = np.sqrt(2) * np.sin(lat_grid)
    
    plt.plot(x_mollweide_grid, y_mollweide_grid, 'k')
    plt.grid(axis='y', color='k')
    
    def add_text(x, lon):
        plt.text(x, 0, str(round(lon * 180 / np.pi, 0)), color='w')
        return
    add_text(np.min(x_mollweide), np.min(lonn))
    add_text(np.min(x_mollweide)/2, np.min(lonn)/2)
    add_text(0, 0)
    add_text(np.max(x_mollweide)/2, np.max(lonn)/2)
    add_text(np.max(x_mollweide), np.max(lonn))
    return

def plot_recv_distribution(full_shadow_vertices, part_shadow_vertices, r_recv_proj):
    """
    输入：
    full_shadow_vertices: 入射孔径内切圆圆心的轨迹在接收板上的投影，构成接收板的全影区域
    part_shadow_vertices: 入射孔径在接收板上的投影，加上以其边缘为圆心的外包络，构成接收板的半影区域
    r_recv_proj: 外包络圆的半径
    输出：
    接收板二维格点的标记阵列
    """
    # 生成二维网格
    x_min, x_max, y_min, y_max = -2, 6, -2, 6
    grid_size = 0.05
    xx, yy = np.meshgrid(np.arange(x_min, x_max, grid_size), np.arange(y_min, y_max, grid_size))
    
    for i in range(len(full_shadow_vertices)):
        full_shadow_vertices[i] = full_shadow_vertices[i][:2]
        part_shadow_vertices[i] = part_shadow_vertices[i][:2]

    # 标记网格点位置关系
    marked_points = np.zeros_like(xx, dtype=int)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            point = (xx[i, j], yy[i, j])
            if  point_polygon_relationship(point,full_shadow_vertices,r=0.01) == 2:
                marked_points[i, j] = 2
            elif point_polygon_relationship(point,part_shadow_vertices,r_recv_proj) == 2 \
                or point_polygon_relationship(point,part_shadow_vertices,r_recv_proj) == 1:
                marked_points[i, j] = 1
            else:
                marked_points[i, j] = 0

    # 绘制标记结果
    cmap = plt.cm.get_cmap('coolwarm', 3)
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(xx, yy, marked_points, cmap=cmap, shading='auto')
    cbar = plt.colorbar()
    cbar.set_ticks([0, 1, 2])
    cbar_labels = ['No Shadow', 'Partial Shadow', 'Full Shadow']
    cbar.set_ticklabels(cbar_labels)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.grid(color='k')

    # 绘制多边形
    polygon_patch_full = plt.Polygon(full_shadow_vertices, closed=True, fill=False, edgecolor='blue', lw=2)
    polygon_patch_part = plt.Polygon(part_shadow_vertices, closed=True, fill=False, edgecolor='green', lw=2)
    plt.gca().add_patch(polygon_patch_full)
    plt.gca().add_patch(polygon_patch_part)
    plt.title('Distribution in Reciving Plane')
    return

estimate_projeffarea(2, [(0,0),(4,0),(2,3)], 1, (2,1), 0.3, 1)