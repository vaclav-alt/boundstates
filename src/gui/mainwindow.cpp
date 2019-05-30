#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "lcp-app.hpp"

#include <iostream>
#include <cstdio>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_btnCalculate_clicked() {
	ui->status->setText("Running");

	LcpApp::Settings s;
	s.a = ui->a->value();
	s.b = ui->b->value();
	s.c = ui->c->value();
	s.eta = ui->eta->value();
	s.ea = ui->EA->value();
	s.mu = ui->mu->value();
	s.E = ui->Emin->value();
	s.NGridPoints = ui->N->value();
	s.basisSize = ui->dvrN->value();
	s.V0file = ui->V0->text().toStdString();
	s.Vlocfile = ui->Vloc->text().toStdString();

	printf("%.10f", s.a);
	LcpApp app(s);
	app.Calculate();

	ui->status->setText("Done");

}
