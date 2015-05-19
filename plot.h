#ifndef _PLOTH_
#define _PLOTH_


#include <stdio.h>
#include <gtk/gtk.h>
#include <glib.h>
#include <glib/gprintf.h>
#include <math.h>

extern struct global_fl_ {
	char zoom_in;
	double x[2];
	double y[2];
} global_fl;

extern struct scale_opt_ {
	int x_range;
	int x_prev;
	int x;
	double y;
} scale_opt;


#define BUFFER_EN_PLACE	2
#define BUFFER_T_PLACE	3
#define BUFFER_dT_PLACE	4

#define ODD(x) ((x % 2) == 0) ? (x=1) : (x=0)
#define OneToZero(x)	(x == 1) ? (x=0) : (x=1)
#define sqr(x)	((x)*(x))

#define IN_EP	0x86
#define OUT_EP	0x02
#define TIMEOUT	10

#define SIZEOF_DATA(x)	((x) ? (1024/4) : (1024))

#define X_NULL	10
#define Y_NULL	20

#define EPS 0.000001
#define CRIT_MAX	0.2


int send_command(const char command, int args);
void calc_results(GtkTextBuffer * txtbuffer, int *a, int *b);

extern char FILETYPE; // 1 for Max spks, 0 for Serg spks 
extern GtkWidget *main_window, *main_statusbar, *coord_statusbar, *entry_range, *entry_delay;

extern int inFile, ok_read;
//char *name_to_read, *name_to_save;

extern pthread_t tid1, tid2;
extern pthread_mutex_t mutex1;
extern int *data[2];
extern double start[2];

extern int rangeset, delay;
extern unsigned char odd_start_draw;

extern cairo_surface_t *surface1;

cairo_surface_t *pt_surface();
double usr_to_win_x(guint width, guint max_x, guint transx, guint x);
double win_to_usr_x(guint width, guint max_x, guint transx, gfloat x);
void plot_bg(cairo_t *cr, int width, int height);
int plot_tics(cairo_t *cr, int width, int height, int *arr, double start_t);
void plot_graph(cairo_t *cr, int width, int height, int *arr);
gboolean graph1_callback (GtkWidget *widget, cairo_t *cr, gpointer user_data);
gboolean graph2_callback (GtkWidget *widget, cairo_t *cr, gpointer user_data);
void *draw(void *user_data);
int max_bubble(int *x, int nums);

#endif
