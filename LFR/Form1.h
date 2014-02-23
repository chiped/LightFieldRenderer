#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/gpu/gpu.hpp>
#include <msclr\marshal_cppstd.h>
#include <math.h>

#using<system.dll>


namespace LFR {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::IO;
	using namespace msclr::interop;
	using namespace cv;
	using namespace std;
	using System::IntPtr;
    using System::Runtime::InteropServices::Marshal;
	
		//Mat Images[16][16];
	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Mat ***Images;
		Mat *Image;

		static double s=0, t=0;
		static int rows, cols;
		static bool quad = true;
		static int gaussianSize = 4;
		static double disparityH = 0, disparityV = 0;
		static int zoom;
	private: System::Windows::Forms::Button^  button2;
	public: 
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::Button^  button4;
	private: System::Windows::Forms::RadioButton^  radioButton1;
	private: System::Windows::Forms::RadioButton^  radioButton2;
	private: System::Windows::Forms::GroupBox^  groupBox1;
	private: System::Windows::Forms::TrackBar^  trackBar1;
	private: System::Windows::Forms::TrackBar^  trackBar2;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TrackBar^  trackBar3;
	private: System::Windows::Forms::Button^  button1;
	public: 
		Form1(void)
		{
			InitializeComponent();
			
			//
			//TODO: Add the constructor code here
			//
			
			System::String^ fileName = "script.txt";
			try 
			{
				Console::WriteLine("trying to open file {0}...", fileName);
				StreamReader^ din = File::OpenText(fileName);

				System::String^ str;
				int count = 0;

				System::String^ delimStr = "#";
				array<Char>^ delimiter = delimStr->ToCharArray( );

				if((str = din->ReadLine()) != nullptr)
					cols = int::Parse(str->Split(delimiter)[1]);
				if((str = din->ReadLine()) != nullptr)
					rows = int::Parse(str->Split(delimiter)[1]);
				
				Images = (Mat***)malloc(rows*sizeof(Mat**));
				for(int k=0;k<rows;k++)
				{
					Images[k] = (Mat**)malloc(cols*sizeof(Mat*));
					for(int l=0;l<cols;l++)
					{
						Images[k][l] = (Mat*)malloc(sizeof(Mat));					
					}
				}
				
				while ((str = din->ReadLine()) != nullptr) 
				{
					count++;
					array<System::String ^>^ splits = str->Split(delimiter);
					for(int i = 0; i < splits->Length; i++)
					this->textBox1->Text = this->textBox1->Text + splits[i] + " ";
					this->textBox1->Text = this->textBox1->Text + "\r\n";

					const char* chars = (const char*)(Marshal::StringToHGlobalAnsi(splits[2])).ToPointer( );
					Images[int::Parse(splits[0])][int::Parse(splits[1])] = new Mat( imread(chars, CV_LOAD_IMAGE_COLOR) );
					Marshal::FreeHGlobal(IntPtr((void*)chars));
				}
				Image = Images[0][0];

				s=t=8;

				//loadImage();
				loadImageDisparity();

			}
			catch (System::Exception^ e)
			{
				if (dynamic_cast<FileNotFoundException^>(e))
					this->textBox1->Text = ("file '{0}' not found " + fileName);
				else
					this->textBox1->Text = ("problem reading file '{0}'" + fileName + e->Message+ " " + e->StackTrace);
			}
			this->ResumeLayout(false);
			this->PerformLayout();
		}

		double getGaussian(int x, int y, double sigma)
		{
			double sigmaSq = sigma * sigma;
			return  1 / (2*CV_PI*sigmaSq) * exp( - ((x-s)*(x-s) + (y-t)*(y-t))/(2*sigmaSq) );
		}

		void gaussianInterpolate()
		{
			int inv = gaussianSize - 1;
			int s0 = (int)floor(s) - inv/2;
			int t0 = (int)floor(t) - inv/2;
			double hD, vD, weight;
			Image = new Mat( Mat::zeros(Images[0][0]->rows, Images[0][0]->cols, Images[0][0]->type()) );
			double total = (gaussianSize - (gaussianSize/4 + 1)); 
			total = total*total*gaussianSize*gaussianSize;
			for(int i = 0; i < gaussianSize; i++)
			{
				for(int j = 0; j < gaussianSize; j++)
				{
					if(s0+inv-i < rows && s0+inv-i >= 0)
						vD = inv - abs(s - (s0 + i));
					else
						vD = inv;
					if(t0+inv-j < cols && t0+inv-j >= 0)
						hD = inv - abs(t - (t0 + j));
					else
						hD = inv;

					if(s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0)
						addWeighted(*Images[s0+i][t0+j], (hD)*(vD)/total, *Image, 1, 0, *Image);
					/*
					if(s0+inv-i < rows && s0+inv-i >= 0)
						weight = getGaussian(s0+i, t0+j, 1);
					else
						weight = getGaussian(s0+i, t0+j, 1) + getGaussian(s0+inv-i, t0+j, 1);

					if (t0+inv-j < cols && t0+inv-j >= 0)
						weight = getGaussian(s0+i, t0+j, 1);
					else
						weight = getGaussian(s0+i, t0+j, 1) + getGaussian(s0+i, t0+inv-j, 1);

					if(s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0)
						addWeighted(*Images[s0+i][t0+j], weight, *Image, 1, 0, *Image);*/
				}
			}
			
		}

		void quadralinearInterpolate()
		{
			int s0 = (int)floor(s);
			int t0 = (int)floor(t);
			double alpha = s - s0;
			double beta = t - t0;

			addWeighted(*Images[s0  ][t0  ], (1-alpha)*(1-beta), 
						*Images[s0  ][t0+1], (1-alpha)*(  beta), 0, *Image);
			addWeighted(*Images[s0+1][t0  ], (  alpha)*(1-beta), *Image, 1, 0, *Image);
			addWeighted(*Images[s0+1][t0+1], (  alpha)*(  beta), *Image, 1, 0, *Image);
/*
			for(int i = 0; i < Image->rows; i++)
				for(int j = 0; j < Image->cols; j++)
					for(int k = 0; k < Image->channels(); k++)
					{

						double sum =      (1-alpha)*(1-beta)*Images[s0  ][t0  ]->at<Vec3b>(i,j)[k]
										+ (1-alpha)*(  beta)*Images[s0  ][t0+1]->at<Vec3b>(i,j)[k]
										+ (  alpha)*(1-beta)*Images[s0+1][t0  ]->at<Vec3b>(i,j)[k]
										+ (  alpha)*(  beta)*Images[s0+1][t0+1]->at<Vec3b>(i,j)[k];

						Image->at<Vec3b>(i,j)[k] = (int)(sum);
					}
					*/
/*
			gpu::GpuMat s0t0;
			s0t0.upload(*Images[s0  ][t0  ]);
			gpu::GpuMat s0t1;
			s0t1.upload(*Images[s0  ][t0+1]);
			gpu::GpuMat s1t0;
			s1t0.upload(*Images[s0+1][t0  ]);
			gpu::GpuMat s1t1;
			s1t1.upload(*Images[s0+1][t0+1]);*/
		}

		void loadImage()
		{
			if(quad)
				quadralinearInterpolate();
			else
				gaussianInterpolate();
			
			IntPtr ptr(Image->ptr());
			Bitmap^ b  = gcnew Bitmap(Image->cols,
								Image->rows,
								Image->step,
								System::Drawing::Imaging::PixelFormat::Format24bppRgb,
								ptr);
			this->pictureBox1->Image = b;

			textBox1->Text = "s = "+s+" t = "+t;
			
		}

		void loadImageDisparity()
		{
			int s0 = (int)floor(s);
			int t0 = (int)floor(t);
			double alpha = s - s0;
			double beta = t - t0;

			Image = new Mat( Mat::zeros(Images[0][0]->rows, Images[0][0]->cols, Images[0][0]->type()) );

			int u0[2][2], v0[2][2];

			for( int u = 0; u < Image->rows; u++)
			{
				for( int v = 0; v < Image->cols; v++)
				{

					for(int i = 0; i < 2; i++)
						for(int j = 0; j < 2; j++)
						{
							u0[i][j] = u + disparityV*(s0+i - s);
							v0[i][j] = v + disparityH*(t0+j - t);
						}
					/*double u1 = u + disparityV*(s0 - s);
					double v1 = v + disparityH*(t0 - t);

					double u2 = u + disparityV*(s0 - s);
					double v2 = v + disparityH*(t0+1 - t);

					double u3 = u + disparityV*(s0+1 - s);
					double v3 = v + disparityH*(t0 - t);

					double u4 = u + disparityV*(s0+1 - s);
					double v4 = v + disparityH*(t0+1 - t);*/

					double value[3] = {0,0,0};

					for(int k = 0; k < 3; k++)
					{
						for(int i = 0; i < 2; i++)
							for(int j = 0; j < 2; j++)
							{
								if(u0[i][j] >= 0 && u0[i][j] <= 239 && v0[i][j] >= 0 && v0[i][j] <= 319)
									value[k] += (1 - abs(s - (s0+i)) )*(1 - abs(t - (t0+j)) ) * Images[s0+i][t0+j]->at<Vec3b>((int)u0[i][j],(int)v0[i][j])[k];
							}

						/*if(u1 >= 0 && u1 <= 239 && v1 >= 0 && v1 <= 319)
							value[k] += (1-alpha)*(1-beta) * Images[s0  ][t0  ]->at<Vec3b>((int)u1,(int)v1)[k];

						if(u2 >= 0 && u2 <= 239 && v2 >= 0 && v2 <= 319)
							value[k] += (1-alpha)*(  beta) * Images[s0  ][t0+1]->at<Vec3b>((int)u2,(int)v2)[k];

						if(u3 >= 0 && u3 <= 239 && v3 >= 0 && v3 <= 319)
							value[k] += (  alpha)*(1-beta) * Images[s0+1][t0  ]->at<Vec3b>((int)u3,(int)v3)[k];

						if(u4 >= 0 && u4 <= 239 && v4 >= 0 && v4 <= 319)
							value[k] += (  alpha)*(  beta) * Images[s0+1][t0+1]->at<Vec3b>((int)u4,(int)v4)[k];*/
					
						Image->at<Vec3b>(u, v)[k] = (int)(value[k]);
					}
				}
			}

			Rect r(zoom, (int)(0.75*zoom), 320 - 2*zoom, (int)(240 - 1.5*zoom));
			Mat temp(*Image, r);
			Mat temp1;
			temp.copyTo(temp1);
			resize(temp1, *Image, Image->size());
			
			IntPtr ptr(Image->ptr());
			Bitmap^ b  = gcnew Bitmap(Image->cols,
								Image->rows,
								Image->step,
								System::Drawing::Imaging::PixelFormat::Format24bppRgb,
								ptr);
			this->pictureBox1->Image = b;

			textBox1->Text = "s = "+s+" t = "+t+"\r\ndisparityH="+disparityH+"\r\ndisparityV="+disparityV;
		}

		void loadImageDisparityGaussian()
		{
			int inv = gaussianSize - 1;
			int s0 = (int)floor(s) - inv/2;
			int t0 = (int)floor(t) - inv/2;
			double hD, vD, weight;
			Image = new Mat( Mat::zeros(Images[0][0]->rows, Images[0][0]->cols, Images[0][0]->type()) );
			double total = (gaussianSize - (gaussianSize/4 + 1)); 
			total = total*total*gaussianSize*gaussianSize;

			double** u0 = new double*[gaussianSize];
			double** v0 = new double*[gaussianSize];
			for(int i = 0; i < gaussianSize; ++i)
			{
				u0[i] = new double[gaussianSize];
				v0[i] = new double[gaussianSize];
			}

			for( int u = 0; u < Image->rows; u++)
			{
				for( int v = 0; v < Image->cols; v++)
				{

					for(int i = 0; i < gaussianSize; i++)
						for(int j = 0; j < gaussianSize; j++)
						{
							u0[i][j] = u + disparityV*(s0+i - s);
							v0[i][j] = v + disparityH*(t0+j - t);
						}
					
					double value[3] = {0,0,0};

					for(int k = 0; k < 3; k++)
					{
						double sum = 0;
						for(int i = 0; i < gaussianSize; i++)
							for(int j = 0; j < gaussianSize; j++)
							{
								if(s0+inv-i < rows && s0+inv-i >= 0)
									vD = inv - abs(s - (s0 + i));
								else
									vD = inv;
								if(t0+inv-j < cols && t0+inv-j >= 0)
									hD = inv - abs(t - (t0 + j));
								else
									hD = inv;
								if(u0[i][j] >= 0 && u0[i][j] <= 239 && v0[i][j] >= 0 && v0[i][j] <= 319 && s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0)
								{
									value[k] += hD*vD * Images[s0+i][t0+j]->at<Vec3b>((int)u0[i][j],(int)v0[i][j])[k];
									sum += hD*vD;
								}

								/*
								
								weight = getGaussian(s0+i, t0+j, 1);

								if(s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0  &&	u0[i][j] >= 0 && u0[i][j] <= 239 && v0[i][j] >= 0 && v0[i][j] <= 319)
								{
									value[k] += weight * Images[s0+i][t0+j]->at<Vec3b>((int)u0[i][j],(int)v0[i][j])[k];
									sum += weight;
								}*/
							}

						Image->at<Vec3b>(u, v)[k] = (int)(value[k]/sum);
					}
				}
			}

			Rect r(zoom, (int)(0.75*zoom), 320 - 2*zoom, (int)(240 - 1.5*zoom));
			Mat temp(*Image, r);
			Mat temp1;
			temp.copyTo(temp1);
			resize(temp1, *Image, Image->size());

			IntPtr ptr(Image->ptr());
			Bitmap^ b  = gcnew Bitmap(Image->cols,
								Image->rows,
								Image->step,
								System::Drawing::Imaging::PixelFormat::Format24bppRgb,
								ptr);
			this->pictureBox1->Image = b;

			textBox1->Text = "s = "+s+" t = "+t;
		}

		void loadImageDisparityZ()
		{
			/*double alpha = s - s0;
			double beta = t - t0;*/

			Image = new Mat( Mat::zeros(Images[0][0]->rows, Images[0][0]->cols, Images[0][0]->type()) );

			int u0[2][2], v0[2][2];

			for( int u = 0; u < Image->rows; u++)
			{
				for( int v = 0; v < Image->cols; v++)
				{
					
					double value[3] = {0,0,0};

					for(int k = 0; k < 3; k++)
					{
						
						double sp, tp,sum=0;

						int up = u;
						int vp = v;
						
						sp = s -0.01*zoom*(up - s);
						tp = t -0.01*zoom*(vp - t);
								
						int s0 = (int)floor(sp);
						int t0 = (int)floor(tp);

						for(int i = 0; i < 2; i++)
							for(int j = 0; j < 2; j++)
							{

								u0[i][j] = up + disparityV*(s0+i - sp);
								v0[i][j] = vp + disparityH*(t0+j - tp);
								double hD = (1 - abs(tp - (t0+j)) );
								double vD = (1 - abs(sp - (s0+i)) );

								if(u0[i][j] >= 0 && u0[i][j] <= 239 && v0[i][j] >= 0 && v0[i][j] <= 319 && s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0)
								{
									value[k] += vD*hD * Images[s0+i][t0+j]->at<Vec3b>((int)u0[i][j],(int)v0[i][j])[k];
									sum += hD*vD;
								}
							}

							/*value[k] += (  alpha)*(  beta) * Images[s0+1][t0+1]->at<Vec3b>((int)u4,(int)v4)[k];*/
					
						Image->at<Vec3b>(u, v)[k] = (int)(value[k]/sum);
					}
				}
			}

			
			IntPtr ptr(Image->ptr());
			Bitmap^ b  = gcnew Bitmap(Image->cols,
								Image->rows,
								Image->step,
								System::Drawing::Imaging::PixelFormat::Format24bppRgb,
								ptr);
			this->pictureBox1->Image = b;

			textBox1->Text = "s = "+s+" t = "+t+"\r\ndisparityH="+disparityH+"\r\ndisparityV="+disparityV;
		}

		void loadImageDisparityGaussianZ()
		{
			int inv = gaussianSize - 1;
			double hD, vD, weight;
			Image = new Mat( Mat::zeros(Images[0][0]->rows, Images[0][0]->cols, Images[0][0]->type()) );
			double total = (gaussianSize - (gaussianSize/4 + 1)); 
			total = total*total*gaussianSize*gaussianSize;

			double** u0 = new double*[gaussianSize];
			double** v0 = new double*[gaussianSize];
			for(int i = 0; i < gaussianSize; ++i)
			{
				u0[i] = new double[gaussianSize];
				v0[i] = new double[gaussianSize];
			}

			for( int u = 0; u < Image->rows; u++)
			{
				for( int v = 0; v < Image->cols; v++)
				{					
					double value[3] = {0,0,0};

					for(int k = 0; k < 3; k++)
					{
						double sum = 0;
						double sp, tp;
						
						sp = s -0.01*zoom*(u - s);
						tp = t -0.01*zoom*(v - t);

						int s0 = (int)floor(sp) - inv/2;
						int t0 = (int)floor(tp) - inv/2;

						for(int i = 0; i < gaussianSize; i++)
							for(int j = 0; j < gaussianSize; j++)
							{
								u0[i][j] = u + disparityV*(s0+i - sp);
								v0[i][j] = v + disparityH*(t0+j - tp);

								vD = inv - abs(sp - (s0+i));
								hD = inv - abs(tp - (t0+j));
								
								if(u0[i][j] >= 0 && u0[i][j] <= 239 && v0[i][j] >= 0 && v0[i][j] <= 319 && s0+i < rows && s0+i >= 0 && t0+j < cols && t0+j >= 0)
								{
									value[k] += hD*vD * Images[s0+i][t0+j]->at<Vec3b>((int)u0[i][j],(int)v0[i][j])[k];
									sum += hD*vD;
								}
							}

						Image->at<Vec3b>(u, v)[k] = (int)(value[k]/sum);
					}
				}
			}

			IntPtr ptr(Image->ptr());
			Bitmap^ b  = gcnew Bitmap(Image->cols,
								Image->rows,
								Image->step,
								System::Drawing::Imaging::PixelFormat::Format24bppRgb,
								ptr);
			this->pictureBox1->Image = b;

			textBox1->Text = "s = "+s+" t = "+t;
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::PictureBox^  pictureBox1;
	private: System::Windows::Forms::TextBox^  textBox1;
	protected: 

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->button4 = (gcnew System::Windows::Forms::Button());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->groupBox1 = (gcnew System::Windows::Forms::GroupBox());
			this->trackBar1 = (gcnew System::Windows::Forms::TrackBar());
			this->trackBar2 = (gcnew System::Windows::Forms::TrackBar());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->trackBar3 = (gcnew System::Windows::Forms::TrackBar());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBox1))->BeginInit();
			this->groupBox1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar2))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar3))->BeginInit();
			this->SuspendLayout();
			// 
			// pictureBox1
			// 
			this->pictureBox1->Location = System::Drawing::Point(63, 42);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(320, 240);
			this->pictureBox1->TabIndex = 0;
			this->pictureBox1->TabStop = false;
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(39, 390);
			this->textBox1->Multiline = true;
			this->textBox1->Name = L"textBox1";
			this->textBox1->ReadOnly = true;
			this->textBox1->ScrollBars = System::Windows::Forms::ScrollBars::Vertical;
			this->textBox1->Size = System::Drawing::Size(361, 130);
			this->textBox1->TabIndex = 1;
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(213, 288);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(26, 23);
			this->button1->TabIndex = 2;
			this->button1->Text = L"▼";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &Form1::button1_Click);
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(389, 151);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(28, 23);
			this->button2->TabIndex = 3;
			this->button2->Text = L"►";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &Form1::button2_Click);
			// 
			// button3
			// 
			this->button3->Location = System::Drawing::Point(33, 151);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(24, 23);
			this->button3->TabIndex = 4;
			this->button3->Text = L"◄";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &Form1::button3_Click);
			// 
			// button4
			// 
			this->button4->Location = System::Drawing::Point(213, 13);
			this->button4->Name = L"button4";
			this->button4->Size = System::Drawing::Size(26, 23);
			this->button4->TabIndex = 5;
			this->button4->Text = L"▲";
			this->button4->UseVisualStyleBackColor = true;
			this->button4->Click += gcnew System::EventHandler(this, &Form1::button4_Click);
			// 
			// radioButton1
			// 
			this->radioButton1->AutoSize = true;
			this->radioButton1->Checked = true;
			this->radioButton1->Location = System::Drawing::Point(6, 21);
			this->radioButton1->Name = L"radioButton1";
			this->radioButton1->Size = System::Drawing::Size(112, 21);
			this->radioButton1->TabIndex = 6;
			this->radioButton1->TabStop = true;
			this->radioButton1->Text = L"Quadralinear";
			this->radioButton1->UseVisualStyleBackColor = true;
			this->radioButton1->CheckedChanged += gcnew System::EventHandler(this, &Form1::radioButton1_CheckedChanged);
			// 
			// radioButton2
			// 
			this->radioButton2->AutoSize = true;
			this->radioButton2->Location = System::Drawing::Point(6, 48);
			this->radioButton2->Name = L"radioButton2";
			this->radioButton2->Size = System::Drawing::Size(89, 21);
			this->radioButton2->TabIndex = 7;
			this->radioButton2->Text = L"Gaussian";
			this->radioButton2->UseVisualStyleBackColor = true;
			this->radioButton2->CheckedChanged += gcnew System::EventHandler(this, &Form1::radioButton2_CheckedChanged);
			// 
			// groupBox1
			// 
			this->groupBox1->Controls->Add(this->radioButton1);
			this->groupBox1->Controls->Add(this->radioButton2);
			this->groupBox1->Location = System::Drawing::Point(425, 369);
			this->groupBox1->Name = L"groupBox1";
			this->groupBox1->Size = System::Drawing::Size(133, 76);
			this->groupBox1->TabIndex = 8;
			this->groupBox1->TabStop = false;
			this->groupBox1->Text = L"Interpolation";
			// 
			// trackBar1
			// 
			this->trackBar1->LargeChange = 1;
			this->trackBar1->Location = System::Drawing::Point(431, 33);
			this->trackBar1->Maximum = -40;
			this->trackBar1->Minimum = -100;
			this->trackBar1->Name = L"trackBar1";
			this->trackBar1->Orientation = System::Windows::Forms::Orientation::Vertical;
			this->trackBar1->Size = System::Drawing::Size(56, 268);
			this->trackBar1->TabIndex = 9;
			this->trackBar1->Value = -40;
			this->trackBar1->Scroll += gcnew System::EventHandler(this, &Form1::trackBar1_Scroll);
			// 
			// trackBar2
			// 
			this->trackBar2->LargeChange = 1;
			this->trackBar2->Location = System::Drawing::Point(487, 33);
			this->trackBar2->Maximum = 7;
			this->trackBar2->Minimum = 1;
			this->trackBar2->Name = L"trackBar2";
			this->trackBar2->Orientation = System::Windows::Forms::Orientation::Vertical;
			this->trackBar2->Size = System::Drawing::Size(56, 268);
			this->trackBar2->TabIndex = 10;
			this->trackBar2->Value = 1;
			this->trackBar2->Scroll += gcnew System::EventHandler(this, &Form1::trackBar2_Scroll);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(36, 369);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(54, 17);
			this->label1->TabIndex = 11;
			this->label1->Text = L"Debug:";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(428, 19);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(46, 17);
			this->label2->TabIndex = 12;
			this->label2->Text = L"Focus";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(480, 19);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(63, 17);
			this->label3->TabIndex = 13;
			this->label3->Text = L"Aperture";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(42, 321);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(48, 17);
			this->label4->TabIndex = 14;
			this->label4->Text = L"Zoom:";
			// 
			// trackBar3
			// 
			this->trackBar3->LargeChange = 1;
			this->trackBar3->Location = System::Drawing::Point(88, 321);
			this->trackBar3->Maximum = 3;
			this->trackBar3->Minimum = -3;
			this->trackBar3->Name = L"trackBar3";
			this->trackBar3->Size = System::Drawing::Size(312, 56);
			this->trackBar3->TabIndex = 15;
			this->trackBar3->Scroll += gcnew System::EventHandler(this, &Form1::trackBar3_Scroll);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(570, 539);
			this->Controls->Add(this->trackBar3);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->trackBar2);
			this->Controls->Add(this->trackBar1);
			this->Controls->Add(this->groupBox1);
			this->Controls->Add(this->button4);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->pictureBox1);
			this->Name = L"Form1";
			this->Text = L"Form1";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBox1))->EndInit();
			this->groupBox1->ResumeLayout(false);
			this->groupBox1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar2))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar3))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
				 s = min(s+0.1, rows-1);
				 //loadImage();
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussian();
			 }
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
			 t = min(t+0.1, cols-1);
			 //loadImage();
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
		 }
private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) {
			 t = max(t-0.1, 0);
			 //loadImage();
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
		 }
private: System::Void button4_Click(System::Object^  sender, System::EventArgs^  e) {
			 s = max(s-0.1, 0);
			 //loadImage();
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
		 }
private: System::Void radioButton1_CheckedChanged(System::Object^  sender, System::EventArgs^  e) {
			 quad = true;
			 //loadImage();
			 loadImageDisparityZ();
		 }
private: System::Void radioButton2_CheckedChanged(System::Object^  sender, System::EventArgs^  e) {
			 quad = false;
			 //loadImage();
			 loadImageDisparityGaussianZ();
		 }
private: System::Void trackBar1_Scroll(System::Object^  sender, System::EventArgs^  e) {
			 disparityH = 0.1* trackBar1->Value;
			 disparityV = 0.1* trackBar1->Value;
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
		 }
private: System::Void trackBar2_Scroll(System::Object^  sender, System::EventArgs^  e) {
			 //gaussianSize = (int)pow(2.0, trackBar2->Value + 1);
			 gaussianSize = 2 + trackBar2->Value*2;
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
			
		 }
private: System::Void trackBar3_Scroll(System::Object^  sender, System::EventArgs^  e) {
			 zoom = trackBar3->Value;
			 if(quad)
				 loadImageDisparityZ();
			 else
				loadImageDisparityGaussianZ();
		 }
};
}

