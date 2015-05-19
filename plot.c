#include "plot.h"


struct global_fl_ global_fl;
struct scale_opt_ scale_opt;

cairo_surface_t *pt_surface()
{
	float red=0.0;
	float green=0.0;
	float blue=0.0;
	
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 6, 6);
  
	cairo_t *context = cairo_create (surface);
	
	cairo_set_line_width (context, 0.005);
	cairo_arc(context, 3, 3, 1.5, 0, 2.0*M_PI);
	cairo_set_source_rgb(context, red, green, blue);
	cairo_stroke_preserve(context);
	cairo_set_source_rgb(context, red, green, blue);
	cairo_fill(context);
	
	return surface;
}

double usr_to_win_x(guint width, guint max_x, guint transx, guint x)
{
	double t = ((double)(width-2.0*X_NULL))/(max_x-1.0-transx)*x + X_NULL - transx*((double)(width-2.0*X_NULL))/(max_x-1.0-transx);

	return t;
}

double win_to_usr_x(guint width, guint max_x, guint transx, gfloat x)
{
	double t = (x-X_NULL)*(max_x-1.0)/(width-2.0*X_NULL) + transx;
	
	return t;
}

void plot_bg(cairo_t *cr, int width, int height)
{
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_set_source_rgba (cr, 0.0, 0.0, 0.0, 1.0);
	cairo_fill(cr);
}

int plot_tics(cairo_t *cr, int width, int height, int *arr, double start_t)
{
	int i=0, max_y = max_bubble(arr, SIZEOF_DATA(FILETYPE)) ;
	double tic_x, mini_tic_x, tic_y;
	char coord_x[4], coord_y[4];
	static guint tics_x_label[] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
	static guint tics_y_label[] = {0, 500, 1000, 1500, 2000, 2050};
	int min_i = (SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 - scale_opt.x;
	int max_i = SIZEOF_DATA(FILETYPE)-1-(SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 + scale_opt.x;
	gfloat x_range = (gfloat) (max_i-min_i-2.0*scale_opt.x);
	
	
	if(min_i < 0 || max_i>SIZEOF_DATA(FILETYPE)-1 || max_i<0) {
		printf("FT = %d, min_i = %d (%d - %d, %d), max_i = %d\n", FILETYPE, min_i, scale_opt.x_range, SIZEOF_DATA(FILETYPE), scale_opt.x, max_i); 
		perror("min i or max i error"); 
		return -1;
	}
	if(max_y<1) max_y = 1;
	
	
	cairo_set_source_rgb (cr, 245.0/255.0, 255.0/255.0, 250.0/255.0);
	cairo_set_line_width (cr, 0.5);
	
	//X axis
	cairo_move_to(cr, 2*X_NULL, height-Y_NULL);
	cairo_line_to(cr, width, height-Y_NULL);
	cairo_stroke(cr);
	
	//Y axis
	cairo_move_to(cr, 2*X_NULL, height-Y_NULL);
	cairo_line_to(cr, 2*X_NULL, Y_NULL);
	cairo_stroke(cr);
	
	for(i=0; i<12; i++) {
		cairo_set_line_width (cr, 2.0);
		tic_x = (double)(2*X_NULL+(width-2.0*X_NULL)*(double)(tics_x_label[i]-min_i)/(x_range));
	//	mini_tic_x = (double)(2*X_NULL+(width-2.0*X_NULL)*(double)(tics_x_label[i]+50-min_i)/(x_range));
		//printf("tic_x = %.2f\n", tic_x);
		
		//if(tic_x < (double)width && mini_tic_x > 0) {
			g_sprintf(coord_x, "%d", tics_x_label[i]);
			cairo_move_to(cr, tic_x, (double)(height-Y_NULL-0.5+4.0));
			cairo_line_to(cr, tic_x, (double)(height-Y_NULL-0.5-4.0));
			cairo_stroke(cr);
		
			cairo_move_to(cr, tic_x-9.0, (double)(height-Y_NULL-0.5+15.0));
			cairo_show_text(cr, coord_x);
			
	//		cairo_set_line_width (cr, 1.0);
	//		cairo_move_to(cr, mini_tic_x, (double)(height-Y_NULL-0.5+3.0));
	//		cairo_line_to(cr, mini_tic_x, (double)(height-Y_NULL-0.5-3.0));
	//		cairo_stroke(cr);
		//	printf("minit tic_x = %.2f\n", tic_x);
		//}
	}
	
	cairo_set_line_width (cr, 2.0);
	for(i=0; i<6; i++) {
		if(tics_y_label[i] >= max_y) break;
		if(i == 5 && fabs(scale_opt.y - 1.0) <= EPS) break;
		tic_y = height-Y_NULL+(2*Y_NULL-height)*(double)(tics_y_label[i]*scale_opt.y/(max_y) - (scale_opt.y-1.0));
		
		if(fabs(tic_y - (height-Y_NULL)) <= EPS) continue;
		//printf("ticy = %.2f, max_y = %.2f\n", tic_y, (double)max_y);
		
		g_sprintf(coord_y, "%d", tics_y_label[i]);
		cairo_move_to(cr, 2*X_NULL+4.0, tic_y);
		cairo_line_to(cr, 2*X_NULL-4.0, tic_y);
		cairo_stroke(cr);
		
		cairo_move_to(cr, 0.0, tic_y-5);
		cairo_show_text(cr, coord_y);
	}
	
	//Start positions
	cairo_set_line_width (cr, .5);
	cairo_move_to(cr, (double)(2*X_NULL+(width-2.0*X_NULL)*(double)(start_t-min_i)/(x_range)), height-Y_NULL);
	cairo_line_to(cr, (double)(2*X_NULL+(width-2.0*X_NULL)*(double)(start_t-min_i)/(x_range)), Y_NULL);
	cairo_stroke(cr);
	
	return 1;
}

void plot_graph(cairo_t *cr, int width, int height, int *arr)
{
	int i, max_y = max_bubble(arr, SIZEOF_DATA(FILETYPE));
	if(max_y<1) max_y = 1;
	
	int min_i = (SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 - scale_opt.x;
	int max_i = SIZEOF_DATA(FILETYPE)-1-(SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 + scale_opt.x;
	gfloat x_range = (gfloat) (max_i-min_i-2.0*scale_opt.x);

	if(min_i < 0 || max_i>SIZEOF_DATA(FILETYPE)-1) return ;
	if(max_y<1) max_y = 1;

	cairo_set_source_rgba (cr, 255.0/255.0, 230.0/255.0, 0.0/255.0, 1.0);
	for(i=min_i; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		cairo_mask_surface(cr, pt_surface(), (double)(2*X_NULL+(width-2.0*X_NULL)*(double)(i-min_i)/(x_range)-3.0), (height-Y_NULL+(2*Y_NULL-height)*(double)( arr[i]*scale_opt.y/(max_y) - (scale_opt.y-1.0))-4.0));
	//	printf("%.2f ", (height-Y_NULL-(2*Y_NULL-height)*(double)arr[i]/max_y)-4.0);
	//	printf("i=%.2f; %d\t", (i-min_i)/x_range, i);
	}
//	printf("ARR[100] = %.2f\n",((double)arr[100]/sca-(double)arr[100]*(1.0/scale_opt.y-1.0))/(max_y));
	
//	printf("\n data[1000] = %d ------------------------------------------\n", arr[1000]);
}

gboolean graph1_callback (GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
	GtkAllocation allocation;
	
	gtk_widget_get_allocation (widget, &allocation);
	
	cairo_set_source_surface (cr, surface1, 0, 0);
	
	plot_bg(cr, allocation.width, allocation.height);
	plot_graph(cr, allocation.width, allocation.height, data[0]);
	
	if(plot_tics(cr, allocation.width, allocation.height, data[0], start[0]) < 0) {
		scale_opt.x = scale_opt.x_prev;
		gtk_widget_queue_draw(main_window);
	}
}

gboolean graph2_callback (GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
	GtkAllocation allocation;
	
	gtk_widget_get_allocation (widget, &allocation);
	
	plot_bg(cr, allocation.width, allocation.height);
	plot_graph(cr, allocation.width, allocation.height, data[1]);
	plot_tics(cr, allocation.width, allocation.height, data[1], start[1]);
}
