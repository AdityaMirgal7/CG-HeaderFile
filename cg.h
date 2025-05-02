//boundaryfill function
void boundaryfill(int x,int y, int fcol,int bcol)
{
	int c=getpixel(x,y);
	if(c!=fcol && c!=bcol)
	{
		putpixel(x,y,fcol);
		boundaryfill(x+1,y,fcol,bcol);
		boundaryfill(x,y+1,fcol,bcol);
		boundaryfill(x-1,y,fcol,bcol);
		boundaryfill(x,y-1,fcol,bcol);
	}	
}


//floodfill function
void floodfill(int x,int y,int oldcol,int newcol)
{
	int c=getpixel(x,y);
	if(c==oldcol)
	{
		putpixel(x,y,newcol);
		floodfill(x+1,y,oldcol,newcol);
		floodfill(x-1,y,oldcol,newcol);
		floodfill(x,y+1,oldcol,newcol);
		floodfill(x,y-1,oldcol,newcol);
	}
}




//dda line function
void dda(int x1, int y1, int x2, int y2)
{

	float dx,dy,i,x,y,step;
	float xinc,yinc;	
	dx=x2-x1;
	dy=y2-y1;
	printf("\n%d\n%d",&dx,&dy);
	

	if(abs(dx)>abs(dy))
	{step=abs(dx);
	}

	else
	{
		step=abs(dy);
	}
xinc=dx/step;
yinc=dy/step;
x=x1;
y=y1;
putpixel(x,y,WHITE);

for(i=0;i<step;i++)
{x=x+xinc;
y=y+yinc;
putpixel(x,y,YELLOW);
delay(50);
}
}


//bresenham line
void bline(int x1, int y1, int x2, int y2)
{
	float dx,dy,x,y,i,e;
	dx=abs(x2-x1);
	dy=abs(y2-y1);
	e=2*(dy-dx);
	x=x1;
	y=y1;
	for(i=0;i<=dx;i++)
	{
	putpixel(x,y,WHITE);
	delay(50);
	
	while(e>=0)
	{y=y+1;
	e=e-2*dx;
	}
	x=x+1;
	e=e+2*dy;
}
}


//bezier curve
int bezier(int x0,int y0,int x1,int y1,int x2,int y2,int x3,int y3)
{

    int i;
	double put_x,put_y;
	double t;
for(t=0.0;t<1.0;t=t+0.001)
	{
		put_x=pow(1-t,3)*x0+pow(1-t,2)*3*t*x1+x2*t*t*(1-t)+x3*pow(t,3);
		put_y=pow(1-t,3)*y0+pow(1-t,2)*3*t*y1+y2*t*t*(1-t)+y3*pow(t,3);
		putpixel(put_x,put_y,WHITE);
		
	}
}

//midpoint line
void drawLine(int x1, int y1, int x2, int y2)
 {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int x = x1, y = y1;
    int incrE = 2 * dy;
    int incrNE = 2 * (dy - dx);
    int d = 2 * dy - dx;

    putpixel(x, y, WHITE);

    while (x < x2) {
        if (d <= 0) {
            d += incrE;
            x++;
        } else {
            d += incrNE;
            x++;
            y++;
        }
        putpixel(x, y, WHITE);
    }  
	  
}

//koch curve

void koch(int x1, int y1, int x2, int y2, int it)
{
 float angle = 60*M_PI/180;
 int x3 = (2*x1+x2)/3;
 int y3 = (2*y1+y2)/3;

 int x4 = (x1+2*x2)/3;
 int y4 = (y1+2*y2)/3;

 int x = x3 + (x4-x3)*cos(angle)+(y4-y3)*sin(angle);
 int y = y3 - (x4-x3)*sin(angle)+(y4-y3)*cos(angle);

 if(it > 0)
 {
  koch(x1, y1, x3, y3, it-1);
  koch(x3, y3, x, y, it-1);
  koch(x, y, x4, y4, it-1);
  koch(x4, y4, x2, y2, it-1);
 }
 else
 {

  line(x1, y1, x3, y3);
  line(x3, y3, x, y);
  line(x, y, x4, y4);
  line(x4, y4, x2, y2);
 }
}

//bresenham circle
void drawCircle(int xc, int yc, int x, int y) {
    putpixel(xc+x, yc+y, YELLOW);
    putpixel(xc-x, yc+y, YELLOW);
    putpixel(xc+x, yc-y, YELLOW);
    putpixel(xc-x, yc-y, YELLOW);
    putpixel(xc+y, yc+x, YELLOW);
    putpixel(xc-y, yc+x, YELLOW);
    putpixel(xc+y, yc-x, YELLOW);
    putpixel(xc-y, yc-x, YELLOW);
}

void bresenhamCircle(int xc, int yc, int r) {
    int x = 0, y = r;
    int d = 3 - 2 * r;

    drawCircle(xc, yc, x, y);
    while (y >= x) {
        x++;

        if (d > 0) {
            y--;
            d = d + 4 * (x - y) + 10;
        } else {
            d = d + 4 * x + 6;
        }

        drawCircle(xc, yc, x, y);
    }
}

//floating point
void drawfCircle(int r)
{
	// Consider a rectangle of size N*N
	int N = 2*r+1;

	int x, y; // Coordinates inside the rectangle

	// Draw a square of size N*N.
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// Start from the left most corner point
			x = i-r;
			y = j-r;

			// If this point is inside the circle, print it
			if (x*x + y*y <= r*r+1 )
				printf(".");
			else // If outside the circle, print space
				printf(" ");
			printf(" ");
		}
		printf("\n");
	}
}


    
    //midpoint ellipse

// Function to draw an ellipse using midpoint algorithm
void drawEllipse(int xc, int yc, int rx, int ry) {
    int x = 0, y = ry;
    int p1 = ry * ry - rx * rx * ry + (rx * rx) / 4;
    
    while (2 * ry * ry * x <= 2 * rx * rx * y) {
        putpixel(xc + x, yc - y, WHITE);
        putpixel(xc - x, yc - y, WHITE);
        putpixel(xc + x, yc + y, WHITE);
        putpixel(xc - x, yc + y, WHITE);
        
        if (p1 < 0) {
            x++;
            p1 = p1 + 2 * ry * ry * x + ry * ry;
        } else {
            x++;
            y--;
            p1 = p1 + 2 * ry * ry * x - 2 * rx * rx * y + ry * ry;
        }
    }

    int p2 = ry * ry * (x + 0.5) * (x + 0.5) + rx * rx * (y - 1) * (y - 1) - rx * rx * ry * ry;
    
    while (y >= 0) {
        putpixel(xc + x, yc - y, WHITE);
        putpixel(xc - x, yc - y, WHITE);
        putpixel(xc + x, yc + y, WHITE);
        putpixel(xc - x, yc + y, WHITE);
        
        if (p2 > 0) {
            y--;
            p2 = p2 - 2 * rx * rx * y + rx * rx;
        } else {
            x++;
            y--;
            p2 = p2 + 2 * ry * ry * x - 2 * rx * rx * y + rx * rx;
        }
    }
}

//neighbours
// Function to find neighbors of a point on a circle using Bresenham algorithm
void findCircleNeighbors(int x_center, int y_center, int radius) {
    int x = radius;
    int y = 0;
    int p = 1 - radius; // Initial decision parameter

    // Print the initial point
    printf("(%d, %d)\n", x_center + x, y_center + y);

    // Print the neighbors using Bresenham algorithm
    while (x > y) {
        y++;

        // Mid-point is inside or on the perimeter of the circle
        if (p <= 0) {
            p = p + 2 * y + 1;
        }
        // Mid-point is outside the perimeter of the circle
        else {
            x--;
            p = p + 2 * y - 2 * x + 1;
        }

        // All the perimeter points have already been printed
        if (x < y) {
            break;
        }

        // Printing the generated point and its reflection in the other octants
        printf("(%d, %d)\n", x_center + x, y_center + y);
        printf("(%d, %d)\n", x_center - x, y_center + y);
        printf("(%d, %d)\n", x_center + x, y_center - y);
        printf("(%d, %d)\n", x_center - x, y_center - y);

        // If the generated point is on the line x = y, the perimeter points have already been printed
        if (x != y) {
            printf("(%d, %d)\n", x_center + y, y_center + x);
            printf("(%d, %d)\n", x_center - y, y_center + x);
            printf("(%d, %d)\n", x_center + y, y_center - x);
            printf("(%d, %d)\n", x_center - y, y_center - x);
        }
    }
}





//cohensutherlanD(Line clipping)
// Define the region codes
#define INSIDE 0x0
#define LEFT 0x1
#define RIGHT 0x2
#define BOTTOM 0x4
#define TOP 0x8

// Define a structure to represent a 2D point
typedef struct {
    int x, y;
} Point;

// Function to compute the region code for a point
int computeRegionCode(Point p, int xmin, int ymin, int xmax, int ymax) {
    int code = INSIDE;
    
    if (p.x < xmin)
        code |= LEFT;
    else if (p.x > xmax)
        code |= RIGHT;
    if (p.y < ymin)
        code |= BOTTOM;
    else if (p.y > ymax)
        code |= TOP;

    return code;
}

// Function to clip a line segment using Cohen-Sutherland algorithm
void cohenSutherlandClip(Point p1, Point p2, int xmin, int ymin, int xmax, int ymax) {
    int code1 = computeRegionCode(p1, xmin, ymin, xmax, ymax);
    int code2 = computeRegionCode(p2, xmin, ymin, xmax, ymax);
    int accept = 0;

    while (1) {
        if ((code1 == 0) && (code2 == 0)) {
            // Both endpoints are inside the window
            accept = 1;
            break;
        } else if (code1 & code2) {
            // Both endpoints are outside the same region
            break;
        } else {
            // Clip the line against the edge
            int codeOut = code1 ? code1 : code2;
            Point newPoint;

            if (codeOut & TOP) {
                newPoint.x = p1.x + (p2.x - p1.x) * (ymax - p1.y) / (p2.y - p1.y);
                newPoint.y = ymax;
            } else if (codeOut & BOTTOM) {
                newPoint.x = p1.x + (p2.x - p1.x) * (ymin - p1.y) / (p2.y - p1.y);
                newPoint.y = ymin;
            } else if (codeOut & RIGHT) {
                newPoint.y = p1.y + (p2.y - p1.y) * (xmax - p1.x) / (p2.x - p1.x);
                newPoint.x = xmax;
            } else if (codeOut & LEFT) {
                newPoint.y = p1.y + (p2.y - p1.y) * (xmin - p1.x) / (p2.x - p1.x);
                newPoint.x = xmin;
            }

            if (codeOut == code1) {
                p1 = newPoint;
                code1 = computeRegionCode(p1, xmin, ymin, xmax, ymax);
            } else {
                p2 = newPoint;
                code2 = computeRegionCode(p2, xmin, ymin, xmax, ymax);
            }
        }
    }

    if (accept) {
        printf("Line (%d, %d) to (%d, %d) is visible after clipping.\n", p1.x, p1.y, p2.x, p2.y);
        // Draw the visible part of the line.
    } else {
        printf("Line is completely outside the window.\n");
        // Line is entirely outside the clipping window.
    }
}


//2d
void twod(int z)
{

			int tx, ty, i, ch, sx, sy, shx, shy, c,c1,s,u;
			float t; 
			int a[4][2],out[4][2];
			int array[4][2] = {{10, 20}, {50, 20}, {50, 60}, {10, 60}};
			
			
			
			line(0, getmaxy() / 2, getmaxx(), getmaxy() / 2);
			line(getmaxx() / 2, 0, getmaxx() / 2, getmaxy());
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=array[u][0];
					a[u][1]=array[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
			while(1){ 			
			printf("\n 1.Translation \n 2.Rotation \n 3.Scaling \n 4.Reflection \n 5.x-shear \n 6.y-shear \n 7.exit\n");
			printf("\n\nEnter ur choice:");
			scanf("%d", &ch);
			setcolor(RED);
			switch (ch)
			{
			case 1: //  Translation
			printf("\nEnter the value of tx and ty:\n");
			scanf("%d%d", &tx, &ty);
			for (i = 0; i < 4; i++)
			{
			out[i][0] = array[i][0] + tx;
			out[i][1] = array[i][1] + ty;
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);

			}
			break;
			
			case 2: //  Rotation in Anticlockwise or Clockwise
			do
			{
				
			printf("\n1.Anticlockwise \n2.Clockwise \n3.Exit");
			printf("Enter ur choice:");
			scanf("%d", &c);
			switch (c)
			{
			case 1: //  Rotation in Anticlockwise
			printf("\nEnter the value of theta:");
			scanf("%f", &t);
			t = (3.14 / 180) * t;
			for (i = 0; i < 4; i++)
			{
			out[i][0] = ((array[i][0]) * cos(t)) - ((array[i][1]) * sin(t));
			out[i][1] = ((array[i][1]) * cos(t)) + ((array[i][0]) * sin(t));
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
				
			}
			break;
			case 2: //  Rotation in Clockwise
			printf("\nEnter the value of theta:");
			scanf("%f", &t);
			t = (3.14 / 180) * t;
			for (i = 0; i < 4; i++)
			{
			out[i][0] = ((array[i][0]) * cos(t)) + ((array[i][1]) * sin(t));
			out[i][1] = ((array[i][1]) * cos(t)) - ((array[i][0]) * sin(t));
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
			break;
			case 3: // For Exit
			break;
			}
			} 
			
			while (c < 3);
			break;
			case 3: // Scaling
			printf("\nEnter the value of sx and sy:");
			scanf("%d%d", &sx, &sy);
			for (i = 0; i < 4; i++)
			{
			
			out[i][0] = array[i][0] * sx;
			out[i][1] = array[i][1] * sy;
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
			break;
			
			case 4:{
				do{
					printf("1. Reflection about x-axis \n2.Reflection about y-axis \n3.Reflection about Origin\n 4.Exit\n");
					scanf("%d",&c1);
					switch(c1){
						case 1:{
							for (i = 0; i < 4; i++)
							{
							out[i][0] = array[i][0];
							out[i][1] = -array[i][1];
							}
							for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
							break;
						}
						case 2:{
							for (i = 0; i < 4; i++)
							{
							out[i][0] = -(array[i][0]);
							out[i][1] = array[i][1];
							}
							for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
							break;
						}
						case 3:{
							for (i = 0; i < 4; i++)
							{
							out[i][0] = -(array[i][0]);
							out[i][1] = -(array[i][1]);
							}
			          		for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
							break;
						}
						case 4:{
							break;
						}
					}
				}while(c1<4);
				break;
			}
			
			case 5: // x-shear
			printf("\nEnter the value of SHX:");
			scanf("%d", &shx);
			for (i = 0; i < 4; i++)
			{
			out[i][0] = (array[i][0]) + (shx * (array[i][1]));
			
			out[i][1] = (array[i][1]);
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
			break;
			
			
			case 6: // y-shear
			printf("Enter the value of SHY:");
			scanf("%d", &shy);
			for (i = 0; i < 4; i++)
			{
			out[i][0] = (array[i][0]);
			out[i][1] = (array[i][1]) + (shy * (array[i][0]));
			}
			for(s=1;s<2;s++){
				for(u=0;u<4;u++){
					a[u][0]=out[u][0];
					a[u][1]=out[u][1];
				}
					line(getmaxx() / 2 + a[0][0], getmaxy() / 2 - a[0][1], getmaxx()/ 2 + a[1][0], getmaxy() / 2 - a[1][1]);
					line(getmaxx() / 2 + a[1][0], getmaxy() / 2 - a[1][1], getmaxx()/ 2 + a[2][0], getmaxy() / 2 - a[2][1]);
					line(getmaxx() / 2 + a[2][0], getmaxy() / 2 - a[2][1], getmaxx()/ 2 + a[3][0], getmaxy() / 2 - a[3][1]);
					line(getmaxx() / 2 + a[3][0], getmaxy() / 2 - a[3][1], getmaxx()/ 2 + a[0][0], getmaxy() / 2 - a[0][1]);
			}
			break;
			
			
			case 7: // exit
			exit(0);
			break;
			
		}
	}
}




//midpoint circle
void midcircle(int x_mid,int y_mid,int radius)
{
int x,y,dp;
x=0;
y=radius;
dp=1-radius;
do
{
putpixel(x_mid+x,y_mid+y,YELLOW);
putpixel(x_mid+y,y_mid+x,YELLOW);
putpixel(x_mid-y,y_mid+x,YELLOW);
putpixel(x_mid-x,y_mid+y,YELLOW);
putpixel(x_mid-x,y_mid-y,YELLOW);
putpixel(x_mid-y,y_mid-x,YELLOW);
putpixel(x_mid+y,y_mid-x,YELLOW);
putpixel(x_mid+x,y_mid-y,YELLOW);
if(dp<0) {
dp+=(2*x)+1;
}
else{
y=y-1;
dp+=(2*x)-(2*y)+1;
}
x=x+1;
}while(y>x);}


//liang barsky
void lb(int x1,int y1,int x2,int y2,int xmin,int ymin,int xmax,int ymax)
{
	int dx,dy,cout,endl;
float	xx1,yy1,xx2,yy2;
float	p[4];
float	q[4];
    float t1,t2, t[4];
    	
dx = x2-x1;
	dy = y2-y1;
	p[0] = -dx;
	p[1] = dx;
	p[2] = -dy;
	p[3] = dy;
	
	
	q[0] = x1-xmin;
	q[1] = xmax-x1;
	q[2] = y1-ymin;
	q[3] = ymax-y1;
	
	for(int i=0; i<4; i++)
	{
		if(p[i] != 0)
		{
			t[i] = (float)q[i]/p[i];
		}
		
		else 
		if(p[i] == 0 && q[i]<0)
			{
					printf("\nLine completly outside the window\n");
			}
			
			
			
		else 
		if(p[i]==0 && q[i] >= 0)
				{
						printf("\nLine is completly inside the window\n");
				
				}
				
						
	}
	if(t[0] > t[2])
	{
		t1 = t[0];
	}
	else{
	   t1 = t[2];
	}
	
	if(t[1] < t[3])
	{
		t2 = t[1];
	}
	else{
	   t2 = t[3];
	}

	
	if(t1 < t2)
	{
		xx1 = x1 + t1*dx;
		xx2 = x1 + t2*dx;
		yy1 = y1 + t1*dy;
		yy2 = y1 + t2*dy;
	
	
		setcolor(WHITE);
		line(xx1,yy1,xx2,yy2);
		
	}
	else{
		printf("\n Line lies outside the window\n");
	}
	
	}
	
//antialised
void drawPixel(int x, int y, float intensity) {
    putpixel(x, y, WHITE * intensity);

}

void drawAntialiasedLine(int x1, int y1, int x2, int y2) 
{
   

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    float intensity;

    if (dx >= dy) {
        for (int x = x1; x <= x2; x++) {
            intensity = (float)(x - x1) / dx;
            int y = y1 + intensity * dy;
            drawPixel(x, y, 1 - intensity);
            drawPixel(x, y + 1, intensity);
        }
    } else {
        for (int y = y1; y <= y2; y++) {
            intensity = (float)(y - y1) / dy;
            int x = x1 + intensity * dx;
            drawPixel(x, y, 1 - intensity);
            drawPixel(x + 1, y, intensity);
        }
    }
}

