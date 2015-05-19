#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <fcntl.h>

#include <cyusb.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include "plot.h"

static cyusb_handle *h = NULL;
int send_command(const char command, int args);
void calc_results(GtkTextBuffer * txtbuffer, int *a, int *b);

char FILETYPE = 1; // 1 for Max spks, 0 for Serg spks 
GtkWidget *main_window, *main_statusbar, *coord_statusbar, *entry_range, *entry_delay;

int inFile = 0, ok_read = 1;
char *name_to_read, *name_to_save;

pthread_t tid1, tid2;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
int *data[2];
double start[2];

int rangeset = 1000, delay = 35;
unsigned char odd_start_draw = 0;

cairo_surface_t *surface1 = NULL;

static struct {
	char *name_to_read;
	int fd_read;
	FILE *fd_Fread;
	int readpos;
	gint sb_context_id;
	
	char *name_to_save;
} files;

int lengthof(char *str)
{
	int i=0;
	
	while(*(str++) != '\0')
		i++;
		
	return ++i;	
}

int f_meander()
{
	static int i = 0;
	static int width = 60, height = 10;
	
	i++;
	if(i/width == 0) return 0;
	else {
		if (i == 2*width-1) i = 0;
		return height;
	}
}

int f_sin()
{
	static int i = 0;
	double T = 10.0;
	i+=5;
	
	return (int)sin(2.0*M_PI*i/T);
}

extern int max_bubble(int *x, int nums)
{
	int i, z;
	
	for(i=1, z = x[0]; i<=nums-1; i++)
		if(x[i] > z)
			z = x[i];
	
	return z;
}

extern int min_bubble_num(int *x, int nums)
{
	int i, z;
	int max_num;
	
	for(i=1, z = x[0]; i<=nums-1; i++)
		if(x[i] < z) {
			z = x[i];
			max_num = i;
		}
	
	return max_num;
}

void swap(int *a, int *b) {
	int x = *b;
	
	*b = *a;
	*a = x;
}


void clicked_on_ga(GtkWidget *widget, GdkEventButton *event, gpointer   data)
{	
	global_fl.zoom_in = 1;
	
	if(event->button == 3) {
		global_fl.zoom_in = 0;
		gtk_widget_queue_draw(main_window);
	}
		
	printf("clicked, zoomin = %d, %d\n", global_fl.zoom_in, (event->state & GDK_BUTTON1_MASK));
}

void unclicked_on_ga(GtkWidget *widget, GdkEventButton *event, gpointer   data)
{
	if((event->state & GDK_BUTTON1_MASK) && (global_fl.zoom_in == 2)) {
		global_fl.zoom_in = 3;
		gtk_widget_queue_draw(main_window);
	}
	
	printf("unclicked, zoomin = %d\n", global_fl.zoom_in);
}

static void motion_in_ga(GtkWidget *widget, GdkEventMotion *event, gpointer   user_data)
{	
	//printf("Motion on GA! (%.2f, %.2f) STATE = %d\n", event->x, event->y, (event->state & GDK_BUTTON2_MASK));
	GtkAllocation allocation;
	
	gtk_widget_get_allocation (widget, &allocation);
	
	int min_i = (SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 - scale_opt.x;
	int max_i = SIZEOF_DATA(FILETYPE)-1-(SIZEOF_DATA(FILETYPE)-scale_opt.x_range)/2 + scale_opt.x;
	gfloat x_range = (gfloat) (max_i-min_i-2.0*scale_opt.x);
	
	int real_i = (int)round( (event->x-2.0*X_NULL)*x_range/(allocation.width-2.0*X_NULL)+min_i );
	
	//printf("real i = %d, x = %.1f\n", real_i, event->x);
	
	if(real_i >= 0)
		gtk_statusbar_push(GTK_STATUSBAR(coord_statusbar), 0, g_strdup_printf("(%d, %d)", real_i, data[GPOINTER_TO_INT(user_data)][real_i]));
}

static void touch_ga(GtkWidget *widget, gpointer   data)
{
	printf("Touch GA!\n");
}

static gboolean graph_configure_event (GtkWidget         *widget,
                          GdkEventConfigure *event,
                          gpointer           data)
{
	printf("%s\n", __FUNCTION__);
	
	GtkAllocation allocation;

  if (surface1)
    cairo_surface_destroy (surface1);

  gtk_widget_get_allocation (widget, &allocation);
  surface1 = gdk_window_create_similar_surface (gtk_widget_get_window (widget),
                                               CAIRO_CONTENT_COLOR,
                                               allocation.width,
                                               allocation.height/2);
}

int control_test(unsigned char bmRequestType) 
{
	int r;
	unsigned char datac[2];
	
	datac[0] = datac[1] = 15;
	
	r = cyusb_control_transfer(h, bmRequestType, 0x00, 0x0000, 0x0000, datac, 2, 0); //wIndex = 0x0000 or LIBUSB_ENDPOINT_IN | IN_EP 0b0000 0000 1000 xxxx
	
	return datac[0];
}

void *write_to_ep(void *user_data)
{	
	static char reset_start = 1;
	
	int x, i, r, transferred = 0;
	unsigned char *buf = (unsigned char *)malloc(2*SIZEOF_DATA(FILETYPE)*sizeof(unsigned char));
	
	if(!access(files.name_to_save, F_OK)) {
		remove(files.name_to_save);
		printf("File was removed\n");
	}
	
	if(!inFile) inFile = open (files.name_to_save, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	
	while(TRUE) {
		//printf("flag=%d\n", ok_read);
		
		if(ok_read) {
			//send_command('r', 0);
			
			for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
					buf[2*i] = 2*f_meander();
					buf[2*i+1] = 4 & 0b00001111;
			}
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 2*SIZEOF_DATA(FILETYPE), &transferred, TIMEOUT);
			
			if(r != 0 || transferred != 2*SIZEOF_DATA(FILETYPE)) {
				if(transferred != 2*SIZEOF_DATA(FILETYPE))
					printf("Error in bulk write transfer, transferred = %d \n", transferred);
				cyusb_error(r);
				cyusb_close();
				exit(0);
			}
			
			
			usleep(20000);
			//odd_start_draw++;
			
			//if(reset_start)
			//	send_command('r', 0);
			OneToZero(reset_start);
			
			r = cyusb_bulk_transfer(h, IN_EP, buf, 2*SIZEOF_DATA(FILETYPE), &transferred, TIMEOUT);
			if(r != 0 || transferred != 2*SIZEOF_DATA(FILETYPE)) {
				if(transferred != 2*SIZEOF_DATA(FILETYPE))
					printf("Error in bulk transfer, transferred = %d \n", transferred);
				cyusb_error(r);
				cyusb_close();
				exit(0);
			}
			printf("trans = %d\n", transferred);
			
			for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
				data[reset_start][i] = buf[2*i]+256*(buf[2*i+1]&0b00001111);
			//	data[reset_start][i] = f_meander();
			//	printf("%u + 256*%u = %d\n", buf[2*i], buf[2*i+1], ((int *)data)[i]);
			}
			//printf("MAX in data = %d\n", max_bubble(data[reset_start], SIZEOF_DATA(FILETYPE)));
			
			for(i=0; i<SIZEOF_DATA(FILETYPE)/2; i++) {
				swap(&(data[reset_start][i]), &(data[reset_start][SIZEOF_DATA(FILETYPE)-i-1]));
			}
			
		//printf("transferred = %d \n-------------------------------------------------------\n", transferred);
		
			if( write(inFile, data[reset_start], sizeof(int)*SIZEOF_DATA(FILETYPE)) == -1 ) {
				perror("Error in read data!");
				exit(1);
			}
		}
	}
	
	free(buf);
	close(inFile);
	
	return 0;
}

static int times_started = 0;
void *read_from_ep(void *user_data)
{
	printf("read form ep\n");
	static char reset_start = 1;
	
	int x, i, j, r, transferred = 0;
	unsigned char *buf = (unsigned char *)malloc(4*SIZEOF_DATA(FILETYPE)*sizeof(unsigned char));
	static const int CONTROL_REQUEST_TYPE_IN = LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_STANDARD | LIBUSB_RECIPIENT_DEVICE;
	
	if(!access(files.name_to_save, F_OK)) {
		remove(files.name_to_save);
		printf("File was removed\n");
	}
	
	if(!inFile) inFile = open (files.name_to_save, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	
	while(TRUE) {
		if(control_test(CONTROL_REQUEST_TYPE_IN) != 1) { continue; }
		if(ok_read) {
			times_started++;
			//OneToZero(reset_start);
			
			memset(buf, 0, 4*SIZEOF_DATA(FILETYPE)*sizeof(unsigned char));
			
			r = cyusb_bulk_transfer(h, IN_EP, buf, 4*SIZEOF_DATA(FILETYPE), &transferred, 0);
		//	printf("trans = %d\n", transferred);
			if(r != 0) {
				if(transferred != 4*SIZEOF_DATA(FILETYPE))
					printf("Error in bulk transfer, transferred = %d \n", transferred);
				cyusb_error(r);
				cyusb_close();
				exit(0);
			}
			
			for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
				data[0][i] = buf[2*i]+256*(buf[2*i+1]&0b00001111);
				data[1][i] = buf[2*i+512]+256*(buf[2*i+1+512]&0b00001111);
			}
			for(i=0; i<SIZEOF_DATA(FILETYPE)/2; i++) {
				swap(&(data[0][i]), &(data[0][SIZEOF_DATA(FILETYPE)-i-1]));
				swap(&(data[1][i]), &(data[1][SIZEOF_DATA(FILETYPE)-i-1]));
			}
			
		//	send_command('r', 0);
		
			if( write(inFile, data[0], sizeof(int)*SIZEOF_DATA(FILETYPE)) == -1 ) {
				perror("Error in read data!");
				printf("fd = %d\n", inFile);
				exit(1);
			}
			
			if( write(inFile, data[1], sizeof(int)*SIZEOF_DATA(FILETYPE)) == -1 ) {
				perror("Error in read data!");
				printf("fd = %d\n", inFile);
				exit(1);
			}
		}
	}
	
	free(buf);
	close(inFile);
}

static void button_read_cb(GtkWidget *widget, gpointer   user_data)
{
	ok_read = 1;
	int err = pthread_create(&tid1, NULL, read_from_ep, NULL);
    if (err != 0)
		printf("\ncan't create clock thread :[%s]", strerror(err));
	
	err = pthread_create(&tid2, NULL, draw, NULL);
    if (err != 0)
		printf("\ncan't create clock thread :[%s]", strerror(err));
}

static void button_write_cb(GtkWidget *widget, gpointer   data)
{
	int err = pthread_create(&tid1, NULL, write_to_ep, NULL);
    if (err != 0)
		printf("\ncan't create clock thread :[%s]", strerror(err));
	
	err = pthread_create(&tid2, NULL, draw, NULL);
    if (err != 0)
		printf("\ncan't create clock thread :[%s]", strerror(err));
}

int send_command(const char command, int args) 
{
	int r, i, transferred = 0;
	unsigned char *buf = (unsigned char *)calloc(256*2, sizeof(unsigned char));
	
	//word #1 = 1; word #2 = porog - установка порога
	//word #1 = 3; word #2 = 1 - сброс (reset)
	//word #1 = 2; word #2 = задержка - установка задержки
	//word #1 = 5; word #2 = 5 - test
	
	switch(command)
	{
		case 's':
			buf[0] = 1;
			buf[2] = args - 256*(args/256);
			buf[3] = args/256;
			printf("rangeset buf[2] = %u buf[3] = %u args = %d\n", buf[2], buf[3], args);
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, TIMEOUT);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in set range!");
				printf("transferred %d bytes\n", transferred);
			}
			free(buf);
			
			
			break;
		case 'r':
			buf[0] = 3;
			buf[2] = 1;
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, 0);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in reset!");
			}
			free(buf);
			
		//	printf("reset\n");
			break;
		case 'w':
			buf[0] = 2;
			buf[2] = args - 256*(args/256);
			buf[3] = args/256;
			
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, TIMEOUT);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in set delay!");
			}
			free(buf);
			
			
			printf("wait buf[2] = %u buf[3] = %u\n", buf[2], buf[3]);
			break;
		case 't':
			buf[0] = 1;
			buf[2] = 5;
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, TIMEOUT);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in test!");
			}
			free(buf);
			
			printf("test");
			break;
		case 'c':
			buf[0] = 5;
			buf[2] = 1;
		//	buf[0] = 1;	моргать светодиодом
		//	buf[2] = 5;
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, TIMEOUT);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in test!");
			}
			free(buf);
			
			printf("coinc on\n");
			break;
		case 'n':
			buf[0] = 5;
			buf[2] = 2;
			
			r = cyusb_bulk_transfer(h, OUT_EP, buf, 256*2, &transferred, TIMEOUT);
			if(r != 0) {
				cyusb_error(r);
				cyusb_close();
				perror("Error in test!");
			}
			free(buf);
			
			printf("coinc off\n");
			break;
		default:
			break;
	}
}


static void button_range_cb(GtkWidget *widget, gpointer   data)
{
	rangeset = atoi(gtk_entry_get_text(GTK_ENTRY(entry_range)));
	
	send_command('s', rangeset);
}

static void button_reset_cb(GtkWidget *widget, gpointer   data)
{
	send_command('r', 0);
}

static void button_wait_cb(GtkWidget *widget, gpointer   data)
{
	delay = atoi(gtk_entry_get_text(GTK_ENTRY(entry_delay)));
	
	send_command('w', delay);
}

static void button_test_cb(GtkWidget *widget, gpointer   data)
{
	send_command('t', 0);
}

static void button_coinc_on_cb(GtkWidget *widget, gpointer   data)
{
	//coincidence
	printf("coincidence\n");
	
	send_command('c', delay);
}

static void button_coinc_off_cb(GtkWidget *widget, gpointer   data)
{
	//none coincidence
	send_command('n', delay);
}

void *draw(void *user_data) 
{
	static int sixty_seconds = 0;
	
	while(TRUE) {
		printf("DRAW! sixty_seconds = %d\n", sixty_seconds);
		sleep(5);
		pthread_mutex_lock( &mutex1 );
		
		
		ok_read = 0;
		
		gtk_widget_queue_draw_area(main_window, 0, 0, 600, 200);
		gtk_widget_queue_draw_area(main_window, 0, 200, 600, 400);
		
	//	sleep(1);
		
		sixty_seconds += 5;
		if(sixty_seconds == 1200) { printf("Exit times_started = %d\n", times_started); pthread_exit(&tid1); pthread_exit(&tid2);}
		
		ok_read = 1;
		
		
		pthread_mutex_unlock( &mutex1 );
	}
}

int draw_data_from_file(GtkTextBuffer *txtbuffer)
{
	gtk_statusbar_push(GTK_STATUSBAR(main_statusbar), files.sb_context_id, g_strdup_printf("File name to read: %s #%d", files.name_to_read, files.readpos));
	
	if(FILETYPE) {
		if( lseek(files.fd_read, SIZEOF_DATA(FILETYPE)*sizeof(int)*files.readpos, SEEK_SET) < 0 ) 
			return -2;
			
		if( read(files.fd_read, data[0], SIZEOF_DATA(FILETYPE)*sizeof(int)) == -1 )
			return -3;
			
		if( read(files.fd_read, data[1], SIZEOF_DATA(FILETYPE)*sizeof(int)) == -1 )
			return -3;
	}
	else {
		int i, j;
		
		if(fseek(files.fd_Fread, ftell(files.fd_Fread), SEEK_SET) != 0)
			return -2;
		
		for(i=0; i<SIZEOF_DATA(FILETYPE); i++) {
			fscanf(files.fd_Fread, "%d %d\n", &j, &data[0][i]);
		}
		for(i=0; i<SIZEOF_DATA(FILETYPE); i++) {
			fscanf(files.fd_Fread, "%d %d\n", &j, &data[1][i]);
		}
		
		int max0 = max_bubble(data[0], SIZEOF_DATA(FILETYPE))
			, max1 = max_bubble(data[1], SIZEOF_DATA(FILETYPE));
		for(i=0; i<SIZEOF_DATA(FILETYPE); i++) {
			data[0][i] = 4000*data[0][i]/max0;
			data[1][i] = 4000*data[1][i]/max0;
		}
	}
		
	gtk_widget_queue_draw(main_window);
	
	calc_results(txtbuffer, data[0], data[1]);
	
	return 0;
}

static void button_readfile_cb(GtkWidget *widget, gpointer user_data)
{
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
	gint i, res;

	dialog = gtk_file_chooser_dialog_new ("Read File",
                                      GTK_WINDOW(main_window),
                                      action,
										"_Cancel",
                                      GTK_RESPONSE_CANCEL,
                                      "_Open",
                                      GTK_RESPONSE_ACCEPT,
                                      NULL);
	chooser = GTK_FILE_CHOOSER (dialog);

	gtk_file_chooser_set_do_overwrite_confirmation (chooser, TRUE);


	res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
	{
		if(files.name_to_read != NULL)
			free(files.name_to_read);
		
		files.name_to_read = (char *)malloc(lengthof(gtk_file_chooser_get_filename (chooser))*sizeof(char));
		
		files.name_to_read = gtk_file_chooser_get_filename (chooser);
		printf("filename choosen = %s %d\n", files.name_to_read, lengthof(gtk_file_chooser_get_filename (chooser)));
		
		if( strncmp(files.name_to_read, "6_1__", strlen("6_1__")) == 0 ) {
			FILETYPE = 0;
		}
		else {
			FILETYPE = 1;
		}
		
		if(FILETYPE) {
			files.fd_read = open(files.name_to_read, O_RDONLY);
		}
		else {
			files.fd_Fread = fopen(files.name_to_read, "r");
		}
		
		files.readpos = 0;
		printf("FT = %d\n", FILETYPE);
		if( draw_data_from_file((GtkTextBuffer *)user_data) != 0 ) {
			if(FILETYPE) {
				close(files.fd_read);
			}
			else {
				fclose(files.fd_Fread);
			}
		}
			
		printf("2\n");
	}

	gtk_widget_destroy (dialog);
}

static void button_readfile_next_cb(GtkWidget *widget, gpointer   user_data)
{
	files.readpos+=2;

	if( draw_data_from_file((GtkTextBuffer *)user_data) != 0 )
			close(files.fd_read);
}

static void button_readfile_prev_cb(GtkWidget *widget, gpointer   user_data)
{
	files.readpos-=2;
	if(files.readpos < 0)
		files.readpos = 0;

	if( draw_data_from_file((GtkTextBuffer *)user_data) != 0 )
			close(files.fd_read);
}

int find_pick_start_stop(int *a, int *max_num, int *min_num)
{
	int i, j;
    int baseline = (a[0] + a[1] + a[2] + a[3] + a[4]) / 5;
     
    for(i=2, j=1; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
        if(a[i] <= 0.95*a[j]) {
		//	printf("i = %d, a[i] = %d, a[j] = %d\n", i, a[i], a[j]);
			if( (a[i+1] <= a[j]) && (a[i+2] <= a[j]) && (a[i] <= baseline) )
                j = i;
                break;
        }
    }
    *max_num = j-7;
    
    for(i=j; i <= SIZEOF_DATA(FILETYPE)-4; i++) {
		if((a[i] >= a[j-7]) &&	(a[i+1]>=a[j-7]) && (a[i]>=baseline)) {
			*min_num = i; 
			return 0;
		}
	}
    
    //*min_num = 150;
    
    return -1;
}

double trap_area(int *a, int max_num, int min_num, double *base)
{
	int i;
	double area = 0.0;
	double baseline = (a[0] + a[1] + a[2] + a[3] + a[4]) / 5.0;
	if(base != NULL) {
		*base = baseline;
	}
	
	for(i=max_num; i<=min_num-1; i++)
		area += 0.5*(a[i]+a[i+1]-2.0*baseline);
		
	return area;
}

// f(x) = ( A*(exp(-((x)-(x0))/a) - exp(-((x)-(x0))/b)) -B*exp(-((x)-(x0))/c) + F )
#define f(x, x0) ( -2.49991703469447e-09*(exp(-((x)-(x0))/115069.384046802) - exp(-((x)-(x0))/1.45088424300257)) - 517.224857987907*exp(-((x)-(x0))/24.9381408089486) + 2506.78894844553 )

//PICK to approx may10_gen_rt=25_dt=500_1sqrt2 #0 (without attenuation) squre = -27702.6
#define NORMAL_SQUARE	(-27702.6)

double search_min_f(double x0, double dx, double start_search)
{
	int i;
	double z = f(start_search, x0), x0min;
	
	printf("x0=%e, dx=%e, start=%e\n", x0, dx, start_search);
	
	for(i=1; i<(int)(6.0/dx); i++) {
		if(f(start_search+i*dx, x0) < z) {
			z = f(start_search+i*dx, x0);
			x0min = start_search + i*dx;
		}
	}
	
	return x0min;
}

double find_start_pick(int *x, char flag)
{
//	#define f(x, x0) ( -0.000725124002756408*(exp(-((x)-(x0))/27461468.6) - exp(-((x)-(x0))/1.31784519694857)) -1855.99996745509*exp(-((x)-(x0))/34.3120374841288) + 3298.00001427814 )
	#define STEP	0.001 
	#define	EPSILON	0.000001
	double x0INIT = 52.0;
	
	double S_optim, S_plus, S_minus, S_test;
	double S_plus_prev, S_minus_prev;
	double x0new = x0INIT;
	int i, j;
	int start_approx = 13;
	int end_approx = 20;
	
	int x_min_num = min_bubble_num(x, SIZEOF_DATA(FILETYPE));
	
	if(flag == 0) {
		double baseline = (x[0] + x[1] + x[2] + x[3] + x[4]) / 5.0;
		int max, min;
		find_pick_start_stop(x, &max, &min);
		
		double normalization = NORMAL_SQUARE/trap_area(x, max, min, NULL);
		printf("normalization = %.4f\n", normalization);
		
		for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
			x[i] = (int)( baseline - normalization*(baseline-(double)x[i]) );
		} 
	}
	
	for(i=0; i<x_min_num; i++) {
		if((x[i] <= 0.97*x[i-1]) && (x[i] <= 0.97*x[i-1])) {
			start_approx = i;
			break;
		}
	}
 
	end_approx = start_approx + 7;
	printf("start = %d; end = %d approx\n", start_approx, end_approx);
	/*
	for(i = x_min_num; i<SIZEOF_DATA(FILETYPE)-1; i++) {
		if( (x[i] >= 0.7*x[x_min_num]) && (x[i+1] >= 0.7*x[x_min_num-4]) ) {
			x0INIT = (double)(i);
			break;
		}
	}
	printf("x_min_num = %d, x0INIT = %.4f\n", x_min_num, x0INIT);
	*/
	S_optim = 0.0;
	S_test = 0.0;
	
	for(i=start_approx; i<=end_approx-1; i++) {
		S_optim += fabs(x[i] - f((double)i, x0INIT));
		S_test += fabs(x[i] - f((double)i, 53.7495));
	}
	
	S_plus = S_minus = 0.0;
	for(j=1; j <= (int)(3.0/STEP); j++) {
		S_plus_prev = S_plus;
		S_minus_prev = S_minus;
		S_plus = S_minus = 0.0;
		
		for(i=start_approx; i<=end_approx-1; i++) {
			S_plus += fabs(x[i] - f((double)i, x0INIT + j*STEP));
			S_minus += fabs(x[i] - f((double)i, x0INIT - j*STEP));
		}
		
		if(S_plus <= S_minus) {
			if(S_plus < S_optim) {
				S_optim = S_plus;
				x0new = x0INIT + STEP*j;
				if(fabs(S_plus - S_plus_prev) <= EPSILON) {
					break;
				}
			}
		}
		else {
			if(S_minus < S_optim) {
		//		printf("j=%d Soptim = %.2f, S+ = %.2f, S- = %.2f, S_test = %.2f; x0new = %.2f\n", j, S_optim, S_plus, S_minus, S_test, x0new);
				S_optim = S_minus;
				x0new = x0INIT - STEP*j;
				if(fabs(S_minus - S_minus_prev) <= EPSILON) {
					break;
				}
			}
		}
	}
	
	printf("j=%d Soptim = %.2f, S+ = %.2f, S- = %.2f, S_test = %.2f; x0new = %.2f\n", j, S_optim, S_plus, S_minus, S_test, x0new);
	
	return x0new;
//	return search_min_f(x0new, STEP/10.0, start_approx);
}


double find_start_cftrace(int *x)
{
	int i;
	for(i=1; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		x[i] = 3000 - x[i];
	}
	
	int *CFTrace = (int *)calloc(SIZEOF_DATA(FILETYPE), sizeof(int));
	for(i=1; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		CFTrace[i] = 0.5*x[i-1] - x[i-1-2];
	}
	FILE *CFTfile = fopen("cftrace.txt", "w");
	for(i=1; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		fprintf(CFTfile, "%d %d\n", i, CFTrace[i]);
	}
	fclose(CFTfile);
	
	double background = (CFTrace[5] + CFTrace[6] + CFTrace[7] + CFTrace[8] + CFTrace[9])/5.0;
	double a, b, x0 = -1.0;
	
	for(i=10; i<SIZEOF_DATA(FILETYPE); i++) {
		if((CFTrace[i-1] >= background) && (CFTrace[i] <= background) && (CFTrace[i-1] >= background+5)) {
			printf("i = %d, background = %.2f, data[5] = %d\n", i, background, CFTrace[5]);
			a = (double)(CFTrace[i]-CFTrace[i-1]);
			b = (double)CFTrace[i] - a*i;
			
			x0 = (background - b)/a;
			break;
		}
	}
	
	free(CFTrace); 

	return x0;
}


/*
double find_start_pick(int *x, char flag)
{
	const int MAX_AMPL = 4000;
	int i, j;
	double z;
	double start_pos;
	
	for(i=1, z = MAX_AMPL-x[0]; i<=SIZEOF_DATA(FILETYPE)-2; i++)
		if (MAX_AMPL-x[i] > z &&  MAX_AMPL-x[i-1] >= 0.2*(MAX_AMPL-x[i]) &&  MAX_AMPL-x[i+1] >= 0.2*(MAX_AMPL-x[i])) {z = MAX_AMPL-x[i]; j = i;}
	z = ((x[j]+x[j+1]+x[j+2])/3.0-(x[0]+x[1]+x[2])/3.0)*CRIT_MAX + (x[0]+x[1]+x[2])/3.0;
	
	if(flag == 0) {
		z = 1300.0;
	}
	else {
		z = 1500.0;
	}
	printf("z = %.2f, j = %d\n", z, j);
	for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++)
		if((double)x[i] <= z) {
			return (start_pos = ( z-(x[i]-(x[i]-x[i-1])*i) ) / (x[i]-x[i-1]) );
		}
}*/

int det3(int a[3][3])
{
	return (a[0][0]*a[1][1]*a[2][2])-(a[0][0]*a[1][2]*a[2][1])
		+(a[0][1]*a[1][2]*a[2][0])-(a[0][1]*a[1][0]*a[2][2])
		+(a[0][2]*a[1][0]*a[2][1])-(a[0][2]*a[1][1]*a[2][0]);
}

double find_start_pick_lsm(int *x)
{
	const int start_n = 160, end_n = 170; // number of points to aproximation
	int i, j;
	int A[3][3], B[3]; // A - matix for Least Square Method
	double X[3];
	double z;
	double start_pos;
	
	for(i=0; i<3; i++) {
		B[i] = 0;
		for(j=0; j<3; j++)
			A[i][j] = 0;
	}
	
	for(i=start_n; i<=end_n; i++) {
		A[0][0] += i*i;
		A[0][1] += i;
		A[1][0] += i*i*i;
		A[2][0] += i*i*i*i;
		
		B[0] += x[i];
		B[1] += i*x[i];
		B[2] += i*i*x[i];
		printf("A00 = %d\n", A[0][0]);
	}

	A[1][1] = A[2][2] = A[0][0];
	A[0][2] = (end_n-start_n+1);
	A[2][1] = A[1][0];
	A[1][2] = A[0][1];

	for(i=0; i<3; i++) {
		printf("B[%d]=%d\n", i, B[i]);
		for(j=0; j<3; j++)
			printf("A[%d][%d]=%d ", i, j, A[i][j]);
		printf("\n");
	}

	int detx[3][3] = {{B[0],A[0][1],A[0][2]},{B[1],A[1][1],A[1][2]},
						{B[2],A[2][1],A[2][2]}};
	int dety[3][3] = {{A[0][0],B[0],A[0][2]},{A[1][0],B[1],A[1][2]},
						{A[2][0],B[2],A[2][2]}};
	int detz[3][3] = {{A[0][0],A[0][1],B[0]},{A[1][0],A[1][1],B[1]},
						{A[2][0],A[2][1],B[2]}};
	
	printf("detA = %d\n", det3(A));
	
	if(det3(A) != 0) {
		X[0] = (double)det3(detx)/det3(A);
		X[1] = (double)det3(dety)/det3(A);
		X[2] = (double)det3(detz)/det3(A);
	}
	for(i=0; i<3; i++)
		printf("X[%d] = %.2f ", i, X[i]);
	printf("\n");
	
	for(i=1, z = x[0]; i<=SIZEOF_DATA(FILETYPE)-2; i++)
		if(x[i] > z && x[i-1] >= 0.2*x[i] && x[i+1] >= 0.2*x[i]) {z = x[i]; j = i;}
	z = ((x[j]+x[j+1]+x[j+2])/3.0-(x[0]+x[1]+x[2])/3.0)*CRIT_MAX + (x[0]+x[1]+x[2])/3.0;
	
	double step = 0.001;
	double xi, xi_next;
	int numOFiter = (end_n-start_n+1)*1000;
	for(i=0; i<=numOFiter; i+=1) {
		xi = i*step+start_n;
		if( ((X[0]*xi*xi + X[1]*xi + X[2] - z) >= EPS)) return xi;
	}
	
	return 0;
}

void calc_results(GtkTextBuffer *txtbuffer, int *a, int *b)
{
	double s_trap[2];
	int max[2], min[2];
	double base[2];
	GtkTextIter iter_start, iter_stop;
	gchar *tempstr;
	
	gtk_text_buffer_get_iter_at_line(txtbuffer, &iter_start, BUFFER_EN_PLACE);
	gtk_text_buffer_get_iter_at_line_offset(txtbuffer, &iter_stop, BUFFER_EN_PLACE+4, 0);
	gtk_text_buffer_delete(txtbuffer, &iter_start, &iter_stop);
	
	find_pick_start_stop(a, &max[0], &min[0]);
//	max[0] -= 5;
//	min[0] = 80;
	s_trap[0] = trap_area(a, max[0], min[0], &base[0]);
	printf("max:%d, min:%d; a[16]=%d\n", max[0], min[0], a[16]);
	
	find_pick_start_stop(b, &max[1], &min[1]);
//	max[1] -= 5;
//	min[1] = 80;
	s_trap[1] = trap_area(b, max[1], min[1], &base[1]);
	printf("strap = %.4f, %.4f\n", s_trap[0], s_trap[1]);
	
	//start[0] = find_start_pick_lsm(a);
	//start[1] = find_start_pick_lsm(b);

	//start[0] = find_start_pick(a, 0);
	//start[1] = find_start_pick(b, 1);
	start[0] = find_start_cftrace(a);
	start[1] = find_start_cftrace(b);
	printf("start = %.2f, %.2f\n", start[0], start[1]);
	
	//find_start_pick_lsm(a);
	
	tempstr = g_strdup_printf("E = %.1fK - %.1f (%d-%d), %.1fK - %.1f (%d-%d)\n", s_trap[0]/1000.0, base[0], max[0], min[0], s_trap[1]/1000.0, base[1], max[1], min[1]);
	gtk_text_buffer_get_iter_at_line(txtbuffer, &iter_start, BUFFER_EN_PLACE);
	gtk_text_buffer_insert_with_tags_by_name(txtbuffer, &iter_start, tempstr, -1, "blue_foreground", NULL);
	g_free(tempstr);
	
	tempstr = g_strdup_printf("T = %.4f, %.4f\n", start[0],  start[1]);							
	gtk_text_buffer_get_iter_at_line(txtbuffer, &iter_start, BUFFER_T_PLACE);
	gtk_text_buffer_insert_with_tags_by_name(txtbuffer, &iter_start, tempstr, -1, "red_foreground", NULL);
	g_free(tempstr);
					
	tempstr = g_strdup_printf("ᐃT = %.4f\n", start[0]-start[1]);		
	gtk_text_buffer_get_iter_at_line(txtbuffer, &iter_start, BUFFER_dT_PLACE);
	gtk_text_buffer_insert_with_tags_by_name(txtbuffer, &iter_start, tempstr, -1, "green_foreground", NULL);
	g_free(tempstr);
}


static void button_read_once_cb(GtkWidget *widget, gpointer   user_data)
{
	int x, i, r, transferred = 0;
	unsigned char *buf = (unsigned char *)malloc(4*SIZEOF_DATA(FILETYPE)*sizeof(unsigned char));

	static const int CONTROL_REQUEST_TYPE_IN = LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_STANDARD | LIBUSB_RECIPIENT_DEVICE;
	static const int CONTROL_REQUEST_TYPE_OUT = LIBUSB_ENDPOINT_OUT | LIBUSB_REQUEST_TYPE_STANDARD | LIBUSB_RECIPIENT_ENDPOINT;
	//control_test(CONTROL_REQUEST_TYPE_OUT);
	control_test(CONTROL_REQUEST_TYPE_IN);
	
/*	r = cyusb_bulk_transfer(h, IN_EP, buf, 4*SIZEOF_DATA(FILETYPE), &transferred, TIMEOUT);
	if(r != 0) {
		cyusb_error(r);
		cyusb_close();
		exit(0);
	}
			
	for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		data[0][i] = buf[2*i]+256*(buf[2*i+1]&0b00001111);
	}
	for(i=0; i<SIZEOF_DATA(FILETYPE)/2; i++) {
		swap((data[0]+i), (data[0]+SIZEOF_DATA(FILETYPE)-i-1));
	}*/
	
	r = cyusb_bulk_transfer(h, IN_EP, buf, 4*SIZEOF_DATA(FILETYPE), &transferred, 0);
		//	printf("trans = %d\n", transferred);
	if(r != 0) {
		if(transferred != 4*SIZEOF_DATA(FILETYPE))
			printf("Error in bulk transfer, transferred = %d \n", transferred);
			cyusb_error(r);
			cyusb_close();
			exit(0);
		}
			
	for(i=0; i<=SIZEOF_DATA(FILETYPE)-1; i++) {
		data[0][i] = buf[2*i]+256*(buf[2*i+1]&0b00001111);
		data[1][i] = buf[2*i+512]+256*(buf[2*i+1+512]&0b00001111);
	}
	for(i=0; i<SIZEOF_DATA(FILETYPE)/2; i++) {
		swap(&(data[0][i]), &(data[0][SIZEOF_DATA(FILETYPE)-i-1]));
		swap(&(data[1][i]), &(data[1][SIZEOF_DATA(FILETYPE)-i-1]));
	}
		
	gtk_widget_queue_draw_area(main_window, 0, 0, 600, 600);
		
	send_command('r', 0);
}

static void scale_zoom_changed(GtkAdjustment *adjust, gpointer       user_data)
{
	//printf("value changed %.2f\n", gtk_adjustment_get_value(adjust));

	scale_opt.x_range = (int) gtk_adjustment_get_value(adjust);
	
	gtk_widget_queue_draw(main_window);
}

static void left_zoom_cb(GtkWidget *widget, gpointer   user_data)
{
	scale_opt.x_prev = scale_opt.x;
	scale_opt.x -= 10;
	
	gtk_widget_queue_draw(main_window);
}

static void right_zoom_cb(GtkWidget *widget, gpointer   user_data)
{
	scale_opt.x_prev = scale_opt.x;
	scale_opt.x += 10;
	
	gtk_widget_queue_draw(main_window);
}

static void yzoom_cb(GtkWidget *widget, gpointer user_data)
{
	scale_opt.y *= 2;
	
	gtk_widget_queue_draw(main_window);
}

static void unzoom_cb(GtkWidget *widget, gpointer   user_data)
{
	scale_opt.x_range = SIZEOF_DATA(FILETYPE);
	scale_opt.x = 0;
	scale_opt.y = 1.0;
	
	gtk_adjustment_set_value(GTK_ADJUSTMENT(user_data), scale_opt.x_range);
	
	gtk_widget_queue_draw(main_window);
}




int main(int   argc, char *argv[])
{
	printf("%d, FT = %d\n", SIZEOF_DATA(FILETYPE), FILETYPE);
	
	int i, r;
	
	for(i=0; i<2; i++) {
		data[i] = (int *)malloc(SIZEOF_DATA(FILETYPE)*sizeof(int *));
		if(data[i] == NULL) {
			printf("Error in calloc data\n");
			exit(1);
		}
	}
	
	scale_opt.x_range = SIZEOF_DATA(FILETYPE);
	scale_opt.x = 0;
	scale_opt.y = 1.0;
	
	files.name_to_save = (char *)malloc(40*sizeof(char));
	if(files.name_to_save == NULL) {
		printf("Error in calloc name_to_save\n");
		exit(1);
	}
	//strcpy(name_to_save, argv[1]);
	files.name_to_save = g_strdup_printf("./spks/%s", argv[1]);
	
	printf("File to save: %s\n", files.name_to_save);
	
	
	GtkWidget *graph_area1, *graph_area2, *main_hbox, *vbox_graph, *button_vbox, *button_read, *button_write, *button_set_range;
	GtkWidget *button_reset,  *button_coinc_on, *button_coinc_off, *button_wait, *button_test, *button_save_to, *button_read_once;
	GtkWidget *table_button, *table_readfile;
	GtkWidget *scale_zoom, *button_rghtzoom, *button_leftzoom, *button_yzoom, *button_unzoom;
	GtkAdjustment *adjust_scale;
	GtkWidget *hr1, *hr2;
	GtkWidget *button_readfile, *button_readfile_next, *image_readfile_next, *button_readfile_prev, *image_readfile_prev;
	GtkWidget *view;
	GtkTextBuffer *buffer_calc;
	
	r = cyusb_open();
	printf("0\n");
	if ( r < 0 ) {
	   perror("Error opening library\n");
	}
	else if ( r == 0 ) {
		perror("No device found!\n");
	}
	printf("1\n");
	if ( r > 1 ) {
		perror("More than 1 devices of interest found. Disconnect unwanted devices\n");
	}
	printf("1.5\n");
	if(r > 0) {
		h = cyusb_gethandle(0);
		if ( cyusb_getvendor(h) != 0x04b4 ) {
			perror("Cypress chipset not detected\n");
			cyusb_close();
		}
		printf("1.55\n");
		r = cyusb_kernel_driver_active(h, 0);
		if ( r != 0 ) {
		   perror("kernel driver active. Exitting\n");
		   cyusb_close();
		}
		printf("1.6\n");
		r = cyusb_claim_interface(h, 0);
		if ( r != 0 ) {
		   perror("Error in claiming interface\n");
		   cyusb_close();
		}
		else printf("Successfully claimed interface\n");
		cyusb_clear_halt(h, OUT_EP);
	}
	printf("2\n");
	
	gtk_init (&argc, &argv);

	main_window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title (GTK_WINDOW (main_window), "Read data from PLIS");
	g_signal_connect (main_window, "destroy", G_CALLBACK (gtk_main_quit), NULL);
	//gtk_window_set_resizable(GTK_WINDOW (main_window), TRUE);

	main_statusbar = gtk_statusbar_new();
	files.sb_context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR (main_statusbar), "Statusbar info");
	gtk_statusbar_push(GTK_STATUSBAR(main_statusbar), files.sb_context_id, g_strdup_printf("File name to save: %s", files.name_to_save));

	coord_statusbar = gtk_statusbar_new();

	main_hbox = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_container_add (GTK_CONTAINER (main_window), main_hbox);


	vbox_graph = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
	gtk_box_pack_start(GTK_BOX(main_hbox), vbox_graph, TRUE, TRUE, 0);
	
	graph_area1 = gtk_drawing_area_new ();
	gtk_widget_add_events(graph_area1, GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK);
	gtk_widget_set_size_request (graph_area1, 600, 200);
	g_signal_connect (graph_area1, "draw",
                    G_CALLBACK (graph1_callback), (gpointer) data);
    g_signal_connect(graph_area1, "button-press-event", 
					G_CALLBACK(clicked_on_ga), NULL);
	g_signal_connect (graph_area1, "motion-notify-event",
					G_CALLBACK (motion_in_ga), GINT_TO_POINTER(0)); 
	g_signal_connect (graph_area1, "button-release-event",
					G_CALLBACK (unclicked_on_ga), NULL);
	g_signal_connect (graph_area1, "configure-event",
					G_CALLBACK (graph_configure_event), NULL); 

	gtk_box_pack_start(GTK_BOX(vbox_graph), graph_area1, TRUE, TRUE, 1);
	
	graph_area2 = gtk_drawing_area_new ();
	gtk_widget_add_events(graph_area2, GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK);
	gtk_widget_set_size_request (graph_area2, 600, 200);
	g_signal_connect (G_OBJECT (graph_area2), "draw",
                    G_CALLBACK (graph2_callback), (gpointer) data);
    /*g_signal_connect(graph_area1, "button-press-event", 
					G_CALLBACK(clicked_on_ga), NULL);*/
	g_signal_connect (graph_area2, "motion-notify-event",
					G_CALLBACK (motion_in_ga), GINT_TO_POINTER(1)); 
	/*g_signal_connect (graph_area1, "button-release-event",
					G_CALLBACK (unclicked_on_ga), NULL);  */

	gtk_box_pack_start(GTK_BOX(vbox_graph), graph_area2, TRUE, TRUE, 1);
	
	gtk_box_pack_start(GTK_BOX(vbox_graph), main_statusbar, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox_graph), coord_statusbar, FALSE, FALSE, 0);
	
	
	button_vbox = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
	gtk_box_pack_start(GTK_BOX(main_hbox), button_vbox, FALSE, FALSE, 0);
	
	
	button_read = gtk_button_new_with_label("Read data");
	g_signal_connect (G_OBJECT (button_read), "clicked",
			G_CALLBACK (button_read_cb), (gpointer) data);
	
	button_write = gtk_button_new_with_label("Write data");
	g_signal_connect (G_OBJECT (button_write), "clicked",
			G_CALLBACK (button_write_cb), (gpointer) NULL);
	
	button_set_range = gtk_button_new_with_label("Установить порог");
	g_signal_connect (G_OBJECT (button_set_range), "clicked",
			G_CALLBACK (button_range_cb), (gpointer) NULL);
	
	button_reset = gtk_button_new_with_label("Reset");
	g_signal_connect (G_OBJECT (button_reset), "clicked",
			G_CALLBACK (button_reset_cb), (gpointer) NULL);
	
	button_test = gtk_button_new_with_label("Тест");
	g_signal_connect (G_OBJECT (button_test), "clicked",
			G_CALLBACK (button_test_cb), (gpointer) NULL);
	
	button_wait = gtk_button_new_with_label("Установить задержку");
	g_signal_connect (G_OBJECT (button_wait), "clicked",
			G_CALLBACK (button_wait_cb), (gpointer) NULL);
	
	button_save_to = gtk_button_new_with_label("Save to file");
/*	g_signal_connect (G_OBJECT(button_save_to), "clicked",
			G_CALLBACK(button_save_to_cb), (gpointer) data);*/

	button_read_once = gtk_button_new_with_label("Read once");
	g_signal_connect (G_OBJECT(button_read_once), "clicked",
			G_CALLBACK(button_read_once_cb), (gpointer) data);
			
	button_coinc_on = gtk_button_new_with_label("Coinc on");
	g_signal_connect (G_OBJECT(button_coinc_on), "clicked",
			G_CALLBACK(button_coinc_on_cb), (gpointer) data);
	
	button_coinc_off = gtk_button_new_with_label("Coinc off");
	g_signal_connect (G_OBJECT(button_coinc_off), "clicked",
			G_CALLBACK(button_coinc_off_cb), (gpointer) data);
	
	entry_range = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(entry_range), g_strdup_printf("%d", rangeset));
	
	entry_delay = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(entry_delay), g_strdup_printf("%d", delay));
	
	adjust_scale = gtk_adjustment_new(SIZEOF_DATA(FILETYPE) + 1.0, 25.0, SIZEOF_DATA(FILETYPE)+1.0, 1, 1.0, 1.0);
	scale_zoom = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adjust_scale);
	gtk_scale_set_value_pos (GTK_SCALE (scale_zoom), GTK_POS_LEFT);
	gtk_scale_set_digits (GTK_SCALE (scale_zoom), 0);
	g_signal_connect (G_OBJECT(adjust_scale), "value-changed",
				G_CALLBACK (scale_zoom_changed), NULL);
	
	button_leftzoom = gtk_button_new_with_label("<");
	g_signal_connect (G_OBJECT(button_leftzoom), "clicked",
			G_CALLBACK(left_zoom_cb), NULL);
			
	button_rghtzoom = gtk_button_new_with_label(">");
	g_signal_connect (G_OBJECT(button_rghtzoom), "clicked",
			G_CALLBACK(right_zoom_cb), NULL);
			
	button_yzoom = gtk_button_new_with_label("+");
	g_signal_connect (G_OBJECT(button_yzoom), "clicked",
			G_CALLBACK(yzoom_cb), NULL);
						
	button_unzoom = gtk_button_new_with_label("Unzoom");
	g_signal_connect (G_OBJECT(button_unzoom), "clicked",
			G_CALLBACK(unzoom_cb), adjust_scale);
	
	hr1 = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
			
	
	table_button = gtk_grid_new();
	gtk_grid_set_column_homogeneous(GTK_GRID(table_button), FALSE);
	gtk_grid_attach(GTK_GRID(table_button), button_read, 0, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_write, 1, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_set_range, 1, 1, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), entry_range, 0, 1, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_reset, 0, 2, 1, 1);
	
	gtk_grid_attach(GTK_GRID(table_button), button_coinc_on, 1, 2, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_coinc_off, 1, 3, 1, 1);
	
	gtk_grid_attach(GTK_GRID(table_button), entry_delay, 0, 3, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_save_to, 0, 4, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_read_once, 1, 4, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), scale_zoom, 0, 5, 2, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_leftzoom, 0, 6, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_rghtzoom, 1, 6, 1, 1);

	gtk_grid_attach(GTK_GRID(table_button), button_yzoom, 2, 6, 1, 1);
	gtk_grid_attach(GTK_GRID(table_button), button_unzoom, 3, 6, 1, 1);
	
	gtk_box_pack_start(GTK_BOX(button_vbox), table_button, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(button_vbox), hr1, FALSE, FALSE, 2);
	
	
	view = gtk_text_view_new ();
	buffer_calc = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));
	gtk_text_buffer_set_text (buffer_calc, "Результаты обработки:\n\nE = \nT = \nᐃT = \n", -1);
	gtk_text_buffer_create_tag (buffer_calc, "blue_foreground",
                              "foreground", "blue", NULL);
    gtk_text_buffer_create_tag (buffer_calc, "red_foreground",
                              "foreground", "red", NULL);
    gtk_text_buffer_create_tag (buffer_calc, "green_foreground",
                              "foreground", "green", NULL);
	
	button_readfile = gtk_button_new_with_label("Read file");
	g_signal_connect (G_OBJECT(button_readfile), "clicked",
			G_CALLBACK(button_readfile_cb), (gpointer) buffer_calc);
	
	button_readfile_next = gtk_button_new();
	image_readfile_next = gtk_image_new_from_file("./go-next-rtl.png");
	gtk_button_set_image(GTK_BUTTON(button_readfile_next), image_readfile_next);
	g_signal_connect (G_OBJECT(button_readfile_next), "clicked",
			G_CALLBACK(button_readfile_next_cb), (gpointer) buffer_calc);
			
	button_readfile_prev = gtk_button_new();
	image_readfile_prev = gtk_image_new_from_file("./go-previous-ltr.png");
	gtk_button_set_image(GTK_BUTTON(button_readfile_prev), image_readfile_prev);
	g_signal_connect (G_OBJECT(button_readfile_prev), "clicked",
			G_CALLBACK(button_readfile_prev_cb), (gpointer) buffer_calc);
	
	
	table_readfile = gtk_grid_new();
	gtk_grid_attach(GTK_GRID(table_readfile), button_readfile, 0, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(table_readfile), button_readfile_next, 2, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(table_readfile), button_readfile_prev, 1, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(table_readfile), view, 0, 1, 4, 2);

	gtk_box_pack_start(GTK_BOX(button_vbox), table_readfile, FALSE, FALSE, 0);

	
	gtk_widget_show_all (main_window);

	gtk_main ();

	return 0;
}






// 1 и 2ой разы считывать
