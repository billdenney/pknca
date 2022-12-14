# dplyr left_join

    Code
      joined <- left_join(myresult, joindf)
    Message <rlang_message>
      Joining, by = "ID"

---

    Code
      joined_manual$result <- dplyr::left_join(joined_manual$result, joindf)
    Message <rlang_message>
      Joining, by = "ID"

---

    Code
      joined <- left_join(myconc, joindf)
    Message <rlang_message>
      Joining, by = "ID"

---

    Code
      joined_manual$data <- left_join(joined_manual$data, joindf)
    Message <rlang_message>
      Joining, by = "ID"

