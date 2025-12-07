# [AC-04] Novel Chemicals Discovery Agent

Accenture
(Komposisi Tim: 2 Machine Learning + 3 Front-End & Back-End)

## 📌 Latar Belakang Masalah

Industri Petrokimia memiliki kebutuhan yang sangat besar untuk menemukan bahan kimia baru yang lebih efisien, ramah lingkungan, dan memiliki nilai ekonomi tinggi. Proses penemuan senyawa kimia baru biasanya membutuhkan waktu lama, biaya besar, serta melibatkan banyak eksperimen laboratorium.

Kehadiran teknologi Artificial Intelligence (AI), khususnya dalam bidang Machine Learning (ML) dan Generative AI, membuka peluang untuk mempercepat proses penemuan produk kimia baru. Dengan memanfaatkan data senyawa kimia yang ada, AI dapat mengusulkan kandidat senyawa atau formula baru yang disesuaikan dengan kriteria produk kimia yang diinginkan.

Use case **Novel Chemicals Discovery Agent** dirancang untuk membantu perusahaan Petrokimia melakukan riset dan eksplorasi produk kimia baru dengan lebih cepat, sistematis, dan berbasis data.

## 🎯 Ruang Lingkup & Tujuan Proyek

### Ruang Lingkup:

1. Membangun sebuah portal aplikasi (web-based) yang dapat menerima input berupa kriteria produk kimia yang diinginkan (misalnya: titik didih, viskositas, kelarutan, stabilitas termal).
2. Menggunakan model AI dengan pendekatan Agentic AI untuk menganalisis kriteria tersebut dan menghasilkan rekomendasi produk kimia potensial, lengkap dengan struktur senyawa, formula, dan deskripsi sifat kimianya.
3. Portal juga dapat menampilkan justifikasi singkat dari hasil rekomendasi.

### Tujuan:

- Memberikan pengalaman kepada Peserta dalam merancang sistem berbasis Agentic AI untuk kasus nyata di industri Petrokimia.
- Mengintegrasikan kemampuan data engineering, data science, dan software development dalam satu proyek.
- Menghasilkan prototype sistem yang dapat menjadi fondasi riset lebih lanjut di bidang chemical informatics.

## 🚀 Hasil yang Diharapkan

1. **Portal aplikasi prototype**, dengan kapabilitas berikut:

   - Dapat menerima input kriteria produk kimia dari user.
   - Menghasilkan daftar rekomendasi produk kimia yang sesuai dengan kriteria.
   - Menampilkan formula, struktur molekul, serta sifat kimia dasar dari rekomendasi tersebut.
   - Menyediakan justifikasi untuk produk kimia yang diusulkan.

2. **Dokumentasi proyek**, yang menjelaskan arsitektur sistem, model AI yang digunakan, dataset yang digunakan, serta limitasi sistem.

## 👥 Komposisi Tim

Proyek ini dikerjakan oleh tim yang terdiri dari:

### 2 Machine Learning Engineers

Ayunda Putri - M251D5X0322`<br>`
Muhammad Fiqri J.A. - M253D5Y1274

### 3 Front-End & Back-End Developers

Itsna Akhdan Fadhil - F012D5Y0871`<br>`
Amrully Arun Hadi - F888D5Y0193`<br>`
Aziz Ridhwan Pratama - F002D5Y0325

## Menjalankan Proyek

Command di bawah mengasumsikan bahwa Anda sedang menggunakan bash-like terminal. (BUKAN COMMAND PROMPT ATAUPUN POWERSHELL DARI WINDOWS!)

### Frontend

```bash
$ cd frontend_js
$ npm install
$ node server.js
```

Frontend akan berjalan di port 3000 secara default.

### Backend

```bash
$ python -m venv .venv
$ source .venv/bin/activate
(.venv) $ python -m pip install --upgrade pip
(.venv) $ python -m pip install -Uvr requirements.txt
(.venv) $ uvicorn src.api_service:app
```

Backend akan berjalan di port 8000 secara default.

## 📚 Sumber Daya Referensi

- [https://medium.com/@aitechtoolbox48/ai-research-teams-just-got-a-superpower-inside-microsofts-discovery-platform-1a90a25fb3b9](https://medium.com/@aitechtoolbox48/ai-research-teams-just-got-a-superpower-inside-microsofts-discovery-platform-1a90a25fb3b9)
- [https://developer.nvidia.com/blog/revolutionizing-ai-driven-material-discovery-using-nvidia-alchemi](https://developer.nvidia.com/blog/revolutionizing-ai-driven-material-discovery-using-nvidia-alchemi)
