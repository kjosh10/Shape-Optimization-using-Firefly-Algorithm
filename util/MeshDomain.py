
class MeshDomain:
    """
    this script is used to generate mesh domain

    Attributes:
        domain: list, a list of domain for the mesh grid
        point_array_Curve: list, a list of points of the curve
        n_x: int, total co-ordinates in x direction
        n_y: int, total co-ordinates in y direction
    """
    def __init__(self,domain, point_array_Curve = [], n_x= 3, n_y=3):
        self.domain = domain
        self.point_array_Curve = point_array_Curve
        self.n_x = n_x
        self.n_y = n_y
        self.points_grid = []
        
    def triMeshInRectArea(self, domain_Mesh):
        """
        method to generate triangular mesh near the Rectangular area
        :param domain_Mesh: list, a list of the domain mesh
        :return pointsTriMeshGlobalDomain: list, a list of the points in the mesh
        """
        #print(domain_Mesh)
        #points_grid = [[(domain_Mesh[0][0]+i/(self.n_vert-1)*(domain_Mesh[1][0]-domain_Mesh[0][0]), domain_Mesh[0][1]+k/(self.n_hor-1)*(domain_Mesh[1][1]-domain_Mesh[0][0]) for i in range(self.n_vert)] for k in range(self.n_hor)]
        
        self.points_grid = [[(domain_Mesh[0][0]+i/(self.n_x-1)*(domain_Mesh[1][0]-domain_Mesh[0][0]), domain_Mesh[0][1]+k/(self.n_y-1)*(domain_Mesh[1][1]-domain_Mesh[0][1])) for i in range(self.n_x)] for k in range(self.n_y)]
        #print(self.points_grid)
        
        #points_Bottom_Layer = points_grid[:self.n_x]
        # print(points_grid)
        pointsTriMeshGlobalDomain = []
        pointsTriMeshLocalDomain = []
        for i in range(self.n_y-1):
            for j in range(self.n_x-1):
                pointsTriMeshLocalDomain.append(self.points_grid[i][j])
                pointsTriMeshLocalDomain.append(self.points_grid[i][j+1])
                pointsTriMeshLocalDomain.append(self.points_grid[i+1][j])
                pointsTriMeshGlobalDomain.append(pointsTriMeshLocalDomain)
                #print(pointsTriMeshLocalDomain)
                pointsTriMeshLocalDomain = []
                pointsTriMeshLocalDomain.append(self.points_grid[i][j+1])
                pointsTriMeshLocalDomain.append(self.points_grid[i+1][j])
                pointsTriMeshLocalDomain.append(self.points_grid[i+1][j+1])
                pointsTriMeshGlobalDomain.append(pointsTriMeshLocalDomain)
                #print(pointsTriMeshLocalDomain)
                pointsTriMeshLocalDomain = []
                
        return pointsTriMeshGlobalDomain
    
    def triMeshInCurve(self, points_Curve):
        """
        method to generate triangular mesh near the Curve area
        :param points_Curve: list, a list of points on the curve
        :param pointsTriMeshGlobalDomain: list, a list of the points in the mesh
        """
        pointsTriMeshGlobalDomain = []
        pointsTriMeshLocalDomain = []
        points_Top_Layer = self.points_grid[0][:self.n_x]
        #print(points_Curve)
        #print(points_Top_Layer)
        #points_grid_Curve_Domain = np.vstack((points_Curve, points_Top_Layer))
        
        for i in range(self.n_x-1):
            
            pointsTriMeshLocalDomain.append(points_Curve[i])
            pointsTriMeshLocalDomain.append(points_Curve[i+1])
            pointsTriMeshLocalDomain.append(points_Top_Layer[i])
            pointsTriMeshGlobalDomain.append(pointsTriMeshLocalDomain)
            #print(pointsTriMeshLocalDomain)
            pointsTriMeshLocalDomain = []
            pointsTriMeshLocalDomain.append(points_Curve[i+1])
            pointsTriMeshLocalDomain.append(points_Top_Layer[i])
            pointsTriMeshLocalDomain.append(points_Top_Layer[i+1])
            pointsTriMeshGlobalDomain.append(pointsTriMeshLocalDomain)
            #print(pointsTriMeshLocalDomain)
            pointsTriMeshLocalDomain = []
            
        return pointsTriMeshGlobalDomain
        
        
    def meshPointsCalc(self):
        """
        method to compute points in the mesh

        :return [cst_Xs, cst_Ys, pts]: list, a list of all the relevant data
        """
        lx = self.domain[1][0]
        if (self.point_array_Curve == [] or self.point_array_Curve == [(0,0),[lx/2,0],[lx,0]]):
            points_CST = self.triMeshInRectArea(self.domain)
            pts = np.array((points_CST))
        else:
            #print("problem")
            domain_rect_area = [(0, pointCalc(self.point_array_Curve, 0.5)[1] + (self.domain[1][1]-pointCalc(self.point_array_Curve, 0.5)[1])/10),(self.domain[1][0], self.domain[1][1])]
            #print(domain_rect_area)
            points_CST_Rect_Area = self.triMeshInRectArea(domain_rect_area)
            points_CST = points_CST_Rect_Area
            #print("Bezier Curve")
            #print(pointCalc(self.point_array_Curve, 0.5)[1])
            t_steps  = np.arange(0,1.000001, (1-0)/(self.n_x-1))
            #print(t_steps)
            points_X_Curve, points_Y_Curve = pointCalc(self.point_array_Curve, t_steps)
            points_Curve = [(points_X_Curve[i], points_Y_Curve[i]) for i in range(points_Y_Curve.shape[0])]
            #print(points_Curve)
            points_CST_Curve = self.triMeshInCurve(points_Curve)
            #(points_CST_Curve)
            points_CST = [points_CST_Curve, points_CST_Rect_Area]
            #print(points_CST.shape)
            pts = np.vstack((points_CST[0], points_CST[1]))
        
        cst_Xs = [(0, pts[i][0][0], pts[i][1][0], pts[i][2][0]) for i in range(pts.shape[0])]
        cst_Ys = [(0, pts[i][0][1], pts[i][1][1], pts[i][2][1]) for i in range(pts.shape[0])]
        return cst_Xs, cst_Ys, pts

