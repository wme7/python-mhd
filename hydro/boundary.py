


class OutflowBoundary:

    def __init__(self):

        self.L_wall = [ [self.x0_wall_1d], [self.x0_wall_2d, self.y0_wall_2d],
                        [self.x0_wall_3d, self.y0_wall_3d, self.z0_wall_3d] ]

        self.R_wall = [ [self.x1_wall_1d], [self.x1_wall_2d, self.y1_wall_2d],
                        [self.x1_wall_3d, self.y1_wall_3d, self.z1_wall_3d] ]


    # --------------------------------------------------------------- 1D
    def x0_wall_1d(self, A, Ng):
        for i in range(Ng):
            A[ i  ] = A[ Ng  ]
    def x1_wall_1d(self, A, Ng):
        for i in range(Ng):
            A[-i-1] = A[-Ng-1]

    # --------------------------------------------------------------- 2D
    def x0_wall_2d(self, A, Ng):
        for i in range(Ng):
            A[ i  ,:] = A[ Ng  ,:]
    def x1_wall_2d(self, A, Ng):
        for i in range(Ng):
            A[-i-1,:] = A[-Ng-1,:]

    def y0_wall_2d(self, A, Ng):
        for i in range(Ng):
            A[:, i  ] = A[:, Ng  ]
    def y1_wall_2d(self, A, Ng):
        for i in range(Ng):
            A[:,-i-1] = A[:,-Ng-1]

    # --------------------------------------------------------------- 3D
    def x0_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[ i  ,:,:] = A[ Ng  ,:,:]
    def x1_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[-i-1,:,:] = A[-Ng-1,:,:]

    def y0_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[:, i  ,:] = A[:, Ng  ,:]
    def y1_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[:,-i-1,:] = A[:,-Ng-1,:]

    def z0_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[:,:, i  ] = A[:,:, Ng  ]
    def z1_wall_3d(self, A, Ng):
        for i in range(Ng):
            A[:,:,-i-1] = A[:,:,-Ng-1]
