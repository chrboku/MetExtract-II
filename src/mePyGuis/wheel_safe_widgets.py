from __future__ import annotations

from PySide6 import QtCore, QtWidgets


class CtrlWheelSpinBox(QtWidgets.QSpinBox):
    """Spin box that only reacts to wheel changes while Ctrl is pressed."""

    def wheelEvent(self, event):
        if event.modifiers() & QtCore.Qt.ControlModifier:
            super().wheelEvent(event)
        else:
            event.ignore()


class CtrlWheelDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    """Double spin box that only reacts to wheel changes while Ctrl is pressed."""

    def wheelEvent(self, event):
        if event.modifiers() & QtCore.Qt.ControlModifier:
            super().wheelEvent(event)
        else:
            event.ignore()


class CtrlWheelComboBox(QtWidgets.QComboBox):
    """Combo box that only reacts to wheel changes while Ctrl is pressed."""

    def wheelEvent(self, event):
        if event.modifiers() & QtCore.Qt.ControlModifier:
            super().wheelEvent(event)
        else:
            event.ignore()


def patch_qt_wheel_sensitive_widgets():
    """Replace Qt widget classes with Ctrl-wheel-safe variants once per process."""

    if getattr(QtWidgets, "_metextract_wheel_widgets_patched", False):
        return

    QtWidgets.QSpinBox = CtrlWheelSpinBox
    QtWidgets.QDoubleSpinBox = CtrlWheelDoubleSpinBox
    QtWidgets.QComboBox = CtrlWheelComboBox
    QtWidgets._metextract_wheel_widgets_patched = True
