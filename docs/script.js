document.querySelectorAll('.dropdown-btn').forEach(button => {
    button.addEventListener('click', function() {
        this.classList.toggle('active');  // Добавляем или убираем класс active для поворота галочки
        let dropdownContent = this.nextElementSibling;
        if (dropdownContent.style.display === "block") {
            dropdownContent.style.display 